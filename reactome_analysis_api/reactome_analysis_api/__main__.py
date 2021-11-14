#!/usr/bin/env python3

import csv
import json
import logging
import uuid
import sys
import os

import connexion
from flask import redirect, request, abort, make_response
from prometheus_client import make_wsgi_app, Counter
from reactome_analysis_api import encoder
from reactome_analysis_utils.reactome_storage import ReactomeStorage, ReactomeStorageException
from reactome_analysis_utils.reactome_logging import get_default_logging_handlers
from werkzeug.middleware.dispatcher import DispatcherMiddleware

app = connexion.App(__name__, specification_dir='./swagger/')
app.app.json_encoder = encoder.JSONEncoder

# show the metrics on a specific path
app_dispatch = DispatcherMiddleware(app, {
    '/metrics': make_wsgi_app()
})

app.add_api('swagger.yaml', arguments={'title': 'REACTOME Analysis Service'})

# log all debug messages to a file
logging.basicConfig(level=logging.DEBUG, handlers=get_default_logging_handlers())

# only log connexion and pika errors
logging.getLogger("connexion").setLevel(logging.ERROR)
logging.getLogger("pika").setLevel(logging.ERROR)
logging.getLogger("rediscluster").setLevel(logging.ERROR)

# log INFO from the util package
logging.getLogger("reactome_analysis_utils").setLevel(logging.INFO)

# set this logger to debug
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)

# count upload errors
UPLOAD_ERRORS = Counter("reactome_api_upload_errors", 
                        "Invalid file uploads", labelnames=["extension"])


def main():
    app.run(port=8080, debug=True)


@app.route("/")
def show_mainpage():
    return redirect("/0.1/ui")

@app.route("/upload", methods=["POST"])
def process_file_upload():
    # test whether the file should be stored or returned
    store_file = request.args.get('store', 'true').lower() == "true"

    # make sure only one file is uploaded
    if len(request.files) != 1:
        abort(400, "Incorrect number of uploaded files. Function requires exactly one file.")

    if "file" not in request.files:
        abort(400, "File must be uploaded as 'file' in the form.")

    # get the uploaded file
    user_file = request.files['file']
    user_filename = user_file.filename

    # initialize the return object
    return_object = {"sample_names": None, "top_identifiers": list(), "n_lines": None}
    return_lines = list()
    n_samples = -1

    # read the file
    try:
        all_lines = [line.decode("UTF-8") for line in user_file.readlines()]
    except Exception as e:
        # get the extension
        (name, extension) = os.path.splitext(user_filename)

        if extension == ".xlsx":
            UPLOAD_ERRORS.labels(extension=extension).inc()
            LOGGER.debug("Excel file upload")
            abort(400, "MS Excel files are not supported. Please save as a text file (txt, csv, or tsv).")
        else:
            LOGGER.info("Invalid file {name} uploaded: {error}".format(name = user_filename, error=str(e)))
            UPLOAD_ERRORS.labels(extension="other").inc()
            abort(400, "Uploaded file is not a text file.")

    # make sure the file is not empty
    if len(all_lines) < 1:
        LOGGER.info("Empty file uploaded.")
        abort(400, "The uploaded file seems to be empty.")

    # guess the delimiter
    delimiter = None
    if "\t" in all_lines[0]:
        delimiter = "\t"
    elif ";" in all_lines[0]:
        delimiter = ";"
    elif "," in all_lines[0]:
        delimiter = ","

    if not delimiter:
        abort(500, "Failed to detect used delimiter")

    try:
        csv_reader = csv.reader(all_lines, delimiter=delimiter)
        header_line = csv_reader.__next__()
        current_line = 1
    except Exception:
        LOGGER.info("Malformatted file encountered.")
        UPLOAD_ERRORS.labels(extension="malformatted csv").inc()
        abort(400, "Malformatted text file. Ensure that quoted fields do not span multiple lines.")

    sample_names = None

    # process each entry
    for line in csv_reader:
        current_line += 1

        if n_samples == -1:
            n_samples = len(line)

            # set the sample names
            if sample_names is None:
                # make sure the sample names are unique
                if len(header_line) == n_samples:
                    sample_names = header_line[1:]
                elif len(header_line) == n_samples - 1:
                    sample_names = header_line
                else:
                    abort(400, "Number of header columns does not match number of samples.")

                # remove empty sample names
                sample_names = [sample_name for sample_name in sample_names if sample_name != ""]
                n_samples = len(sample_names)

                # make sure the sample names are unique
                if len(sample_names) != len(set(sample_names)):
                    abort(400, "Duplicate sample names detected. All sample names (labels in the first line) "
                                "must be unique")

                # save the sample names
                return_object["sample_names"] = sample_names

            # make sure the file was parsed more or less correctly
            if n_samples < 2:
                abort(400, "Failed to parse the file. Only one column detected.")

            # start creating the converted object by adding the header line
            return_lines.append("\t" + "\t".join(sample_names))

        # make sure the number of samples is OK
        if len(line) - 1 < n_samples:
            abort(400, "Different number of entries in line {}. File contains {} columns but line {} contains {}"
                  .format(str(current_line), str(n_samples), str(current_line), str(len(line))))

        # make sure the line (= gene / protein id) is not empty
        if len(line[0].strip()) == 0:
            continue

        # save the first few identifiers as samples
        if current_line < 10:
            return_object["top_identifiers"].append(line[0])

        # save the line
        return_lines.append("\t".join(line[0:n_samples+1]))

    # save the results
    return_object["n_lines"] = current_line

    # create the complete result string
    result_string = "\n".join(return_lines)

    # add the file if it shouldn't be saved
    if not store_file:
        return_object["data"] = result_string
    else:
        # store the file
        try:
            storage = ReactomeStorage()

            # create an identifier
            token = "rqu_" + str(uuid.uuid1())

            while storage.request_token_exists(token):
                token = "rqu_" + str(uuid.uuid1())

            # save the data - expire after 6 hours
            storage.set_request_data(token=token, data=result_string, expire=60*60*6)

            return_object["data_token"] = token
        except ReactomeStorageException as e:
            LOGGER.error("Failed to store request data: " + str(e))
            abort(500, "Failed to store request data. Please try again later.")

    # return the JSON data
    response_object = make_response(json.dumps(return_object))
    # Using the content-type "text/html" instead of the more
    # appropriate "application/json" to circumvent the lacking
    # support for JSON in GWT (used by Reactome's pathway browser)
    response_object.headers["Content-Type"] = "text/html"

    return response_object


if __name__ == '__main__':
    main()
