"""
Class that is listening to report generation requests and
creates the actual report file
"""

import json
import logging
import multiprocessing
import os
import sys
import smtplib
from email.headerregistry import Address
from email.message import EmailMessage

import prometheus_client
import rpy2.rinterface as ri
import rpy2.rinterface_lib
import rpy2.robjects as ro
from reactome_analysis_api.models.report_status import ReportStatus
from reactome_analysis_api.models.report_status_reports import ReportStatusReports
from reactome_analysis_utils import reactome_mq, reactome_storage
from reactome_analysis_utils.models import report_request


# create the logger
LOGGER = logging.getLogger(__name__)

# initialize R
ri.initr()

# disable R messages
def ignore_message(message: str) -> None:
    pass

def log_r_warning(message: str) -> None:
    LOGGER.warn("R Error: " + message)

# exit the process if R needs any input on the console
def exit_process(message: str) -> None:
    LOGGER.error("R Console Read triggered: " + message)
    sys.exit(1)

rpy2.rinterface_lib.callbacks.consolewrite_print = ignore_message
rpy2.rinterface_lib.callbacks.consolewrite_warnerror = log_r_warning
rpy2.rinterface_lib.callbacks.consoleread = exit_process

# set the counters
RUNNING_REPORTS = prometheus_client.Gauge("reactome_report_running",
                                          "Number of reports currently being created.")
COMPLETED_REPORTS = prometheus_client.Counter("reactome_report_completed",
                                              "Number of successfully completed reports.")
SENT_EMAILS = prometheus_client.Counter("reactome_report_emails_sent",
                                        "Number of sent e-mails.")


class ReactomeAnalysisReportGenerator:
    def __init__(self):
        """
        Initializes basic member variables
        """
        self._mq = None
        self._storage = None
        self._performed_analyses = 0

    def _get_mq(self):
        """
        Return the current connection to the message queue
        :return: A ReactomeMQ object
        """
        if not self._mq:
            try:
                self._mq = reactome_mq.ReactomeMQ(queue_name=reactome_mq.REPORT_QUEUE)
            except Exception as e:
                LOGGER.error("Failed to connect to MQ service: " + str(e))
                raise Exception("Failed to connect to MQ service.", e)

        return self._mq

    def _get_storage(self):
        """
        Returns the current connection to the reactome storage
        :return: A ReactomeStorage object
        """
        if not self._storage:
            try:
                self._storage = reactome_storage.ReactomeStorage()
            except Exception as e:
                LOGGER.error("Failed to connect to storage service: " + str(e))
                raise Exception("Failed to connect to storage service", e)

        return self._storage

    def shutdown(self):
        """
        Gracefully shutdown all services
        """
        LOGGER.debug("Shutting down gracefully.")
        if self._mq:
            LOGGER.debug("Closing MQ connection.")
            self._mq.close()

    def start_listening(self):
        """
        Connects to the report queue and starts listening for new report creation
        requests. This is a blocking function that will only exit on failure.
        """
        mq = self._get_mq()

        LOGGER.debug("Listening for messages...")
        mq.process_analysis(self._on_new_report)

    def _acknowledge_message(self, channel, method):
        """
        Acknowledges a message and thereby marks the report as complete.
        :param channel: The channel on which the message was received
        :param method: The method object passed to the analysis
        :return:
        """
        # This function is only here to increase the readability of the code.
        channel.basic_ack(delivery_tag=method.delivery_tag)

        # Decrement the number of running reports
        RUNNING_REPORTS.dec()

    def _set_status(self, analysis_id, status: str, completion: float = 0,
                    description: str = None, reports: list = None):
        """
        Set the status for the report
        :param analysis_id: The report's id
        :param status: Current status ("running", "complete", or "failed")
        :param completion: Relative amount of completion
        :param description: Text describing the current status
        :param reports: A list of created reports
        """
        self._get_storage().set_status(
            analysis_identifier=analysis_id, data_type="report",
            status=json.dumps(ReportStatus(
                id=analysis_id, status=status, completed=completion, description=description, reports=reports)
                              .to_dict()))

    def _on_new_report(self, ch, method, properties, body):
        """
        Triggered by a new request for the report generation
        :param ch: The channel the message was received on
        :param method: Method details
        :param properties: Message properties
        :param body: The actual message body (= JSON encoded analysis request)
        """
        LOGGER.debug("Received message.")

        RUNNING_REPORTS.inc()

        # create the request object
        try:
            LOGGER.debug("Decoding JSON string")
            # decode the JSON information
            request = report_request.from_json(body)
        except Exception as e:
            # This means that the application has a major problem - this should never happen
            LOGGER.critical("Failed to create report request object: " + str(e))
            LOGGER.debug("Error details:", exc_info=1)
            # remove the message from the queue
            self._acknowledge_message(ch, method)

            return

        # get the result data
        LOGGER.debug("Creating report for " + request.analysis_id)
        storage = self._get_storage()
        analysis_result = storage.get_result(analysis_identifier=request.analysis_id)

        if analysis_result is None:
            LOGGER.debug("Failed to retrieve result for " + request.analysis_id)
            self._set_status(analysis_id=request.analysis_id, status="failed",
                             description="Failed to retrieve analysis result")
            self._acknowledge_message(ch, method)
            return

        # indicate that the report is being created
        self._set_status(analysis_id=request.analysis_id, status="running", description="Creating reports",
                         completion=0.1)

        # launch the process to create the report
        on_complete_event = multiprocessing.Event()
        report_result_queue = multiprocessing.Queue()

        report_process = ReportGenerationProcess(analysis_result=analysis_result, report_request=request,
                                                 on_complete=on_complete_event, result_queue=report_result_queue)

        LOGGER.debug("Starting report generation process...")
        report_process.start()

        # wait for the process to be completed
        while report_process.is_alive() and not on_complete_event.is_set():
            # test whether the analysis should be interrupted
            if self._get_mq().get_is_shutdown():
                LOGGER.debug("Shutdown triggered, terminating analysis process")
                report_process.terminate()
                report_process.join(0.1)
                return

            # wait for 1 sec - the process will take a while
            self._get_mq().sleep(1)

        LOGGER.debug("Report process finished. on_complete_event = " + str(on_complete_event.is_set()))

        # update the status
        self._set_status(analysis_id=request.analysis_id, status="running", description="Fetching reports",
                         completion=0.8)

        # for potential cleanup
        report_process.join(1)

        # get all results
        generated_results = list()
        # the base_url is taken of an environmental parameter
        base_url = os.getenv("BASE_URL", "http://193.62.55.4")

        while report_result_queue.qsize() > 0:
            report_filename = report_result_queue.get(block=True, timeout=0.5)

            if isinstance(report_filename, Exception):
                # indicate the something went wrong
                LOGGER.error("Report generation exception: " + str(report_filename))
                self._set_status(analysis_id=request.analysis_id, status="running",
                                 description="Failed to create one report", completion=0.9)
            else:
                # store the result file
                if os.path.isfile(report_filename):
                    extension = os.path.splitext(report_filename)[1]

                    if extension == ".xlsx":
                        data_type = "report"
                    elif extension == ".pdf":
                        data_type = "pdf_report"
                    elif extension == ".r":
                        data_type = "r_script"
                    else:
                        LOGGER.error("Unknown extension encountered: " + extension)
                        continue

                    # load the file and save it in storage
                    with open(report_filename, "rb") as binary_reader:
                        LOGGER.debug("Saving " + extension + " report...")
                        binary_data = binary_reader.read()
                        storage.set_result(analysis_identifier=request.analysis_id, result=binary_data,
                                           data_type=data_type)

                        # save the generated report as a ReportStatusReports object
                        generated_results.append(ReportStatusReports(
                            name=ReactomeAnalysisReportGenerator.get_name_for_extension(extension),
                            url="{base_url}/0.1/result/{analysis_id}{extension}".format(
                                base_url=base_url, analysis_id=request.analysis_id, extension=extension),
                            mimetype=ReactomeAnalysisReportGenerator.get_mimetype_for_extension(extension))
                        )

                    # delete the file
                    os.remove(report_filename)

        # send the e-mail if set
        if request.user_mail:
            # get the available visualizations
            available_visualizations = dict()
            analysis_result_obj = json.loads(analysis_result)
            if "reactome_links" in analysis_result_obj and analysis_result_obj["reactome_links"]:
                available_visualizations = \
                    dict([(link["name"], link["url"]) for link in analysis_result_obj["reactome_links"]])

            LOGGER.debug("Sending report e-mail...")
            ReactomeAnalysisReportGenerator.send_email(user_address=request.user_mail,
                                                       available_reports=generated_results,
                                                       available_visualizations=available_visualizations)

        self._acknowledge_message(ch, method)

        # update the status
        if len(generated_results) > 0:
            self._set_status(analysis_id=request.analysis_id, status="complete", completion=1,
                             description="Report generation complete.", reports=generated_results)
            COMPLETED_REPORTS.inc()
        else:
            self._set_status(analysis_id=request.analysis_id, status="failed", completion=1,
                             description="Report generation failed")

        LOGGER.debug("Report creation completed.")

    @staticmethod
    def send_email(user_address, available_reports, available_visualizations: dict):
        """
        Send an e-mail notifying the user about the new
        reports of the analysis.
        :param user_address: The user's e-mail address
        :param available_reports: A list of the extensions of the available reports
        :param available_visualizations: List containing the names (key) and links (value) of available Reactome
                                         visualizations
        """
        # only send an e-mail if at least one report was created
        if len(available_reports) < 1:
            LOGGER.error("Not sending an e-mail as no reports were generated")
            return

        # only try to send an e-mail if the configuration is present
        smtp_server = os.getenv("SMTP_SERVER", None)
        smtp_port = os.getenv("SMTP_PORT", None)

        if not smtp_server or not smtp_port:
            return

        # Create the base text message.
        msg = EmailMessage()
        msg['Subject'] = "Reactome Analysis Complete"
        msg['From'] = Address("Reactome Analysis Service", addr_spec=os.getenv("FROM_ADDRESS", "no-reply@reactome.org"))

        try:
            msg['To'] = user_address
        except Exception:
            LOGGER.info(f"Invalid mail address ${user_address}. Not sending reports.")
            return

        # Create the plain-text content
        plain_report_lines = ["  * {vis_name} (Reactome visualization): {vis_link}".format(vis_name=vis_name,
                                                                  vis_link=available_visualizations[vis_name])
                              for vis_name in available_visualizations]
        plain_report_lines += ["  * {report_name}: {url}".format(
            report_name=report.name,
            url=report.url)
            for report in available_reports]

        msg.set_content("""\
        Dear user,
        
        Your analysis request to the Reactome Gene Set Analysis Service is complete. 
        
        You can download your results here:
        
          {reports}

        Kind regards,
        The Reactome Team

        ----
        DISCLAIMER

        This message contains confidential data/information and is intended 
        only for the individual named. If you are not the named addressee, 
        you should not disseminate, distribute or copy this email. Please 
        notify the sender immediately by email if you have received this email 
        by mistake and delete this email from your system. Email transmission 
        cannot be guaranteed to be secure or error-free, as information could 
        be intercepted, corrupted, lost, destroyed, arrive late or incomplete, 
        or contain viruses. The sender, therefore, does not accept liability for 
        any errors or omissions in the contents of this message which arise as a 
        result of email transmission.
        """.format(reports="\n".join(plain_report_lines)))

        # Add the html version.  This converts the message into a multipart/alternative
        # container, with the original text message as the first part and the new html
        # message as the second part.
        # Create the HTML content
        html_report_lines = ["  <li><a href=\"{vis_link}\">{vis_name} (Reactome visualization)</a></li>"
                             .format(vis_name=vis_name,
                                     vis_link=available_visualizations[vis_name])
                             for vis_name in available_visualizations]
        html_report_lines += ["<li><a href=\"{url}\">{report_name}</a></li>"
            .format(
                url=report.url, report_name=report.name) for report in available_reports]

        msg.add_alternative("""\
        <html>
          <head></head>
          <body>
            <p>Dear user,</p>
            <p>Your analysis request to the Reactome Gene Set Analysis Service 
               is complete. You can download your results here:</p>
            <ul>
              {reports}
            </ul>
            <p>
              Kind regards, <br />
              The Reactome Team
            </p>
            <p>
            <b>Disclaimer:</b></br>
            This message contains confidential data/information and is intended 
            only for the individual named. If you are not the named addressee, 
            you should not disseminate, distribute or copy this email. Please 
            notify the sender immediately by email if you have received this email 
            by mistake and delete this email from your system. Email transmission 
            cannot be guaranteed to be secure or error-free, as information could 
            be intercepted, corrupted, lost, destroyed, arrive late or incomplete, 
            or contain viruses. The sender, therefore, does not accept liability for 
            any errors or omissions in the contents of this message which arise as a 
            result of email transmission.
            </p>
          </body>
        </html>
        """.format(reports="\n".join(html_report_lines)), subtype='html')

        # get the user and password for the smtp server
        smtp_userfile = os.getenv("EMAIL_USER_FILE", None)
        if not smtp_userfile:
            return

        with open(smtp_userfile, "r") as reader:
            smtp_user = reader.read().strip()

        smtp_pwdfile = os.getenv("EMAIL_PASSWORD_FILE", None)
        if not smtp_pwdfile:
            return

        with open(smtp_pwdfile, "r") as reader:
            smtp_pwd = reader.read().strip()

        # Send the message
        try:
            with smtplib.SMTP(host=smtp_server, port=smtp_port) as s:
                s.ehlo()
                s.starttls()
                s.ehlo()
                s.login(user=smtp_user, password=smtp_pwd)
                s.ehlo()
                s.send_message(msg)

            SENT_EMAILS.inc()
        except Exception as e:
            LOGGER.error("Failed to send e-mail to {address}".format(address=user_address))

    @staticmethod
    def get_name_for_extension(extension: str) -> str:
        """
        Returns the report name for the passed extension
        :param extension: The extension
        :return: The matching name
        """
        if extension == ".pdf":
            return "PDF Report"
        elif extension == ".xlsx":
            return "MS Excel Report (xlsx)"
        elif extension == ".r":
            return "R Script"
        else:
            return "Report"

    @staticmethod
    def get_mimetype_for_extension(extension: str) -> str:
        """
        Returns the mimetype for the passed extension
        :param extension: The extension
        :return: The matching mimetype
        """
        if extension == ".pdf":
            return "application/pdf"
        elif extension == ".xlsx":
            return "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        elif extension == ".r":
            return "text/plain"
        else:
            return None


class ReportGenerationProcess(multiprocessing.Process):
    """
    Class used to create the (R-based) report
    in a separate process.
    """
    def __init__(self, analysis_result: str, report_request: report_request.ReportRequest, on_complete: multiprocessing.Event,
                 result_queue: multiprocessing.Queue):
        """
        Initializes a ReportGenerationProcess
        :param analysis_result: The analysis result as a JSON-encoded string
        :param report_request: The report request object
        :param on_complete: Even triggered once the analysis is complete.
        :param result_queue: Queue that will receive the filenames of the result files.
        """
        super().__init__()

        self.analysis_result = analysis_result
        self.report_request = report_request
        self.on_complete = on_complete
        self.result_queue = result_queue

    def run(self) -> None:
        try:
            # inject the analysis_result into the R session
            ri.globalenv["analysis_result_json"] = ri.StrSexpVector([self.analysis_result.decode()])

            # inject the metadata
            ri.globalenv["include_interactors"] = ri.BoolSexpVector([self.report_request.include_interactors])
            ri.globalenv["include_disease"] = ri.BoolSexpVector([self.report_request.include_disease])

            # create the analysis result object
            LOGGER.debug("Creating result R object...")
            ro.reval("""
                library(ReactomeGSA)

                # convert the reactome object
                result_obj <- jsonlite::fromJSON(analysis_result_json)
                reactome_obj <- ReactomeGSA:::convert_reactome_result(result_obj)
            """)

            # create the Excel file
            LOGGER.debug("Creating Excel file ...")

            excel_filename = "/tmp/result_" + self.report_request.analysis_id + ".xlsx"
            self.create_excel_file(excel_filename)

            self.result_queue.put(excel_filename)

            # create the PDF report
            LOGGER.debug("Creating PDF report...")

            pdf_filename = "/tmp/result_" + self.report_request.analysis_id + ".pdf"
            self.create_pdf_report(pdf_filename)

            self.result_queue.put(pdf_filename)

            # create the R script
            LOGGER.debug("Creating R script...")

            r_filename = "/tmp/result_" + self.report_request.analysis_id + ".r"
            self.create_r_script(r_filename)

            self.result_queue.put(r_filename)
        except Exception as e:
            # put the error message in the queue
            LOGGER.error("Error during report generation: " + str(e))
            self.result_queue.put(e)

        finally:
            LOGGER.debug("Setting on_complete")
            self.on_complete.set()

    def create_excel_file(self, filename: str) -> None:
        """
        Create the Excel result file based on the `reactome_obj` in the R session
        :param filename: Path to the Excel file that will be created
        """

        # inject the result path
        ri.globalenv["excel_result_file"] = ri.StrSexpVector([filename])

        ro.reval("""
            # get the pathways table
            pathway_result <- pathways(reactome_obj)

            # create the workbook
            library(openxlsx)

            wb <- createWorkbook(creator = "ReactomeGSA", title = "ReactomeGSA Analysis Result", subject = "Pathway analysis")

            # bold headers
            boldHeader <- createStyle(textDecoration = 'bold')

            # write the pathways
            addWorksheet(wb, 'Pathways')
            writeData(wb, 'Pathways', pathway_result, headerStyle = boldHeader)
            setColWidths(wb, 'Pathways', cols = 1:ncol(pathway_result), widths = 'auto')

            # store all previous sheet names to ensure that they are unique
            previous_sheet_names <- c()

            # add the expression values for every result
            if ("fold_changes" %in% result_types(reactome_obj)) {
                for (dataset_name in names(reactome_obj)) {
                    fold_changes <- get_result(reactome_obj, type = "fold_changes", name = dataset_name)

                    # add the fold-changes to the Excel file
                    sheet_name <- paste0(dataset_name, " - fold changes")

                    # make the sheet name save
                    sheet_name <- gsub("[^A-Za-z0-9_.]", "_", sheet_name)
                    sheet_name <- gsub("_+", "_", sheet_name)

                    if (nchar(sheet_name) > 25) {
                        sheet_name <- paste0(substr(sheet_name, 1, 25), "...")
                    }

                    # ensure that the sheetName is unique
                    counter <- 0
                    base_sheet_name <- sheet_name

                    while (sheet_name %in% previous_sheet_names) {
                        counter <- counter + 1
                        sheet_name <- paste0(base_sheet_name, "-", counter)
                    }

                    previous_sheet_names <- c(previous_sheet_names, sheet_name)

                    # add the sheet to the workbook
                    addWorksheet(wb, sheet_name)
                    writeData(wb, sheet_name, fold_changes, headerStyle = boldHeader)

                    # nice column widths
                    setColWidths(wb, sheet_name, cols = 1:ncol(fold_changes), widths = 'auto')
                }
            }

            # save the workbook
            saveWorkbook(wb, excel_result_file, overwrite = TRUE)
        """)

    def create_pdf_report(self, pdf_filename) -> None:
        """"
        Create the PDF report based on the `reactome_obj` in the
        current R session.
        :param pdf_filename: The target filename of the PDF report. Will be overwritten if it exists.
        """
        # inject the path used to store the xlsx file
        ri.globalenv["pdf_result_file"] = ri.StrSexpVector([pdf_filename])

        try:
            ro.reval("""
                library(ReactomeGSA.report)
                
                options(tinytex.verbose = TRUE)
    
                # create the report
                create_pdf_report(reactome_obj, pdf_result_file, include_disease = include_disease, include_interactors = include_interactors)
            """)
        except Exception as e:
            LOGGER.error("RRuntimeError: " + str(e))

    def create_r_script(self, r_filename: str) -> None:
        """
        Create an R script to load the analysis result into an R session.
        :param r_filename: The target filename of the script
        """
        # the base_url is taken of an environmental parameter
        base_url = os.getenv("BASE_URL", "http://193.62.55.4")

        with open(r_filename, "w") as writer:
            writer.write("""
# This script downloads your recent ReactomeGSA result
# into an R session
#
# Note: The result is only stored for a certain period of time
#       on the ReactomeGSA servers. Therefore, it is highly
#       recommended to store the result locally.

# install the ReactomeGSA package if not available
if (!require(ReactomeGSA)) {{
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

    BiocManager::install("ReactomeGSA")
}}

# load the package
library(ReactomeGSA)

# load the analysis result
result <- get_reactome_analysis_result(analysis_id = "{analysis_id}", reactome_url = "{base_url}")

# save the result
saveRDS(result, file = "my_ReactomeGSA_result.rds")

# get the overview over all pathways
all_pathways <- pathways(result)
""".format(analysis_id=self.report_request.analysis_id, base_url=base_url))
