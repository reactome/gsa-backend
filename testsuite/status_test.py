import logging
import requests
import json
import time
import click
from submit_request import send_request


logger = logging.getLogger(__name__)


@click.command(help="Continuously send status requests to the specific ReactomeGSA service")
@click.option("--url", help="URL of the Reactome service to use")
@click.argument("id")
def test_status(url: str, id: str):
    logging.basicConfig(level=logging.DEBUG)
    urllib_logger = logging.getLogger("urllib3")
    urllib_logger.setLevel(logging.ERROR)

    while True:
        request_url = url + "/0.1/status/" + id
        request_obj = requests.get(request_url)

        if request_obj.status_code != 200:
            logger.error("Error code " + str(request_obj.status_code))
            continue

        data = json.loads(request_obj.content)
        logger.info(data["description"])
        time.sleep(1)


if __name__ == "__main__":
    test_status()
