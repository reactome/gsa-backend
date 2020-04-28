# ReactomeGSA Test Suite

This is an automated test suite for the ReactomeGSA analysis system. It centers around
the `submit_request.py` python script. This script can send arbitrary requests to
a ReactomeGSA analysis system and run automatic tests on the results.

## Usage

To submit a request simply call:

```bash
python3 submit_request.py requests/melanoma_b_cell_example_request.json
```

## Specifying ReactomeGSA systems

The script automatically loads available ReactomeGSA installations from the `sources.json`
file. The file must be located in the same directory as the `submit_request.py` file.

The `sources.json` file contains a (JSON-formatted) list of ReactomeGSA installations:

```json
[
    {
        "name": "reactome-prod",
        "url": "https://gsa.reactome.org"
    },
    {
        "name": "local",
        "url": "http://my-installation.local:12345"
    }
]
```

The second entry shows how you can specify custom ports for the ReactomeGSA installation to test.

If no additional command line parameters are passed, the script will display a prompt to
enter the name of the installation to use. Alternatively, the name can be supplied
using the `-s` / `--server` command line parameter (this parameter must supply the name).

## Request file format

Request specification files are JSON formatted text files that
contain a request object as defined in the 
[ReactomeGSA API specification](https://gsa.reactome.org).

Additionally, these requests objects may contain a "tests" property containing
the specification of tests that are automatically run:

```json
Test Specification:
   'tests': [
       {
           'name': Test name, # name of the test to run
           'dataset': str, # name of the dataset to evaluate
           'type': ["pathways", "fold_changes"] # whether to test number of pathways or fold changes
           'value': numeric # The expected value
       }
   ]
```