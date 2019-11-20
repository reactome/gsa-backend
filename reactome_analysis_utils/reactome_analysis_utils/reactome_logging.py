"""
Collection of logging functions used in the Reactome GSA
service.

At the moment, this only contains a method to generate an
e-mail handler based on the existing e-mail configuration
"""

import logging
import logging.handlers
import os, sys


LOGGER = logging.getLogger(__name__)


class ReactomeSMTPHandler(logging.handlers.MemoryHandler):
    """
    SMTPHandler that sends one single e-mail containing all entries
    in the buffer.
    """
    def __init__(self, capacity=50, flushLevel=logging.ERROR):
        logging.handlers.MemoryHandler.__init__(self, capacity, flushLevel)

        # indicates whether all mail settings were found
        self.has_mail_support = False

        # get the e-mail settings from the environment variables
        user_file = os.getenv("EMAIL_USER_FILE")
        password_file = os.getenv("EMAIL_PASSWORD_FILE")

        if not user_file or not password_file:
            LOGGER.debug("Missing e-mail configuration to create e-mail logger.")
            return

        user = ReactomeSMTPHandler._load_secret_from_file(user_file)
        password = ReactomeSMTPHandler._load_secret_from_file(password_file)

        if not user or not password:
            LOGGER.debug("Missing e-mail credentials to create e-mail logger.")
            return

        # create the handler
        self.smtp_server =  os.environ["SMTP_SERVER"]
        self.smtp_port = os.environ["SMTP_PORT"]
        self.smtp_from = os.environ["FROM_ADDRESS"]
        self.smtp_to = "jgriss@ebi.ac.uk"
        self.smtp_user = user
        self.smtp_password = password

        # set a sensibel formatter
        self.setFormatter(logging.Formatter("[%(asctime)s] %(levelname)s [%(name)s.%(funcName)s:%(lineno)d] %(message)s"))

    def flush(self):
        """
        Flush all messages to a single e-mail
        """
        if len(self.buffer) > 0:
            try:
                # create the nicely formatted string
                messages = list()

                for log_msg in self.buffer:
                    messages.append(self.format(log_msg))

                # create the e-mail message
                from email.message import EmailMessage

                msg = EmailMessage()
                msg['Subject'] = "ReactomeGSA - Error"
                msg['From'] = self.smtp_from
                msg['To'] = self.smtp_to

                msg.set_content("\n".join(messages))

                # send the mail
                import smtplib

                with smtplib.SMTP(host=self.smtp_server, port=self.smtp_port) as s:
                    s.ehlo()
                    s.starttls()
                    s.ehlo()
                    s.login(user=self.smtp_user, password=self.smtp_password)
                    s.ehlo()
                    s.send_message(msg)
            except Exception as e:
                LOGGER.debug("Failed to send log mail: " + str(e))

            self.buffer = []

    @staticmethod
    def _load_secret_from_file(file):
        """
        Load the contents of the file into a string and return it. Returns None if the file
        does not exist.
        :param file: Path to the file. If None None is returned
        :return: The contents (stripped) as string or None if the file does not exist
        """
        if not file or not os.path.isfile(file):
            return None

        with open(file, "r") as reader:
            contents = reader.read()

            return contents.strip()