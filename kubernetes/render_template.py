"""render_template.py

Renders the Kubernetes YAML specification for the REACTOME Analysis
System.

Usage:
    render_template.py --template TEMPLATE_FILE --config CONFIG_FILE --output RESULT_FILE [--k_1_18]
    render_template.py (--help | --usage)

Parameters:
    -t, --template TEMPLATE_FILE            The template file to use. This must be a jinja2
                                            formatted template.
    -c, --config CONFIG_FILE                A YAML formatted file containing all parameters
                                            required for the template
    -o, --output OUTPUT_FILE                Path to the output file (must not exist).
    --k_1_18                                Adapt the template to be compatible with k8s >= 
                                            version 1.18
    -h, --help                              Displays this help
"""

from jinja2 import Environment,FileSystemLoader
import yaml
import base64
import string
import random
from docopt import docopt
import os
import sys
import tempfile
import shutil


def random_password(length=30, debug=False, alphanumeric=False):
    """
    Generate a random string first
    """
    if debug:
        password = "test"
    else:
        letters = string.ascii_lowercase + string.ascii_uppercase + string.digits
        if not alphanumeric:
            letters +=  "/&*+i§![]=-"
        password = ''.join(random.choice(letters) for i in range(length))

    base64_password = base64.encodebytes(password.encode()).decode()

    return base64_password


def fix_1_18(filename: str) -> None:
    """Fix the k8s yaml file to be compatible with k8s version >= 1.18

    :param filename: The file to fix
    :type filename: str
    """
    # create a temporary file
    output_file = tempfile.mktemp(suffix=".yaml")

    with open(filename, "r") as reader:
        with open(output_file, "w") as writer:
            for line in reader:
                line = line.replace("extensions/v1beta1", "apps/v1")

                writer.write(line)

    # replace the original file
    shutil.move(output_file, filename)


def main():
    args = docopt(__doc__)

    template_file = args["--template"]
    config_file = args["--config"]
    output_file = args["--output"]
    update_1_18 = args["--k_1_18"]

    # make sure all files exist
    for input_file in [template_file, config_file]:
        if not os.path.isfile(input_file):
            print("Error: {} does not exist".format(input_file))
            sys.exit(1)
    # make sure the output file doesn't
    if os.path.isfile(output_file):
        print("Error: {} already exists".format(output_file))
        sys.exit(1)

    # set the environment and add the password function
    env = Environment(loader=FileSystemLoader(os.path.dirname(template_file)))
    env.globals.update(random_password=random_password)

    # load the config
    with open(config_file, "r") as config_reader:
        conf = yaml.safe_load(config_reader)

    # make sure the namepsace has a valid value
    if not "namespace" in conf or len(conf["namespace"].strip()) == 0:
        conf["namespace"] = None

    if not "storage_class" in conf or len(conf["storage_class"].strip()) == 0:
        conf["storage_class"] = None

    if not "proxy" in conf["reactome_worker"] or len(conf["reactome_worker"]["proxy"].strip()) == 0:
        conf["reactome_worker"]["proxy"] = None

    if not "debug" in conf:
        conf["debug"] = False

    if not "mail_error_address" in conf:
        conf["mail_error_address"] = ""

    t = env.get_template(os.path.basename(template_file))
    t.stream(conf).dump(output_file)

    # update to k8s 1.18 if set
    if update_1_18:
        fix_1_18(output_file)


if __name__ == "__main__":
    main()

