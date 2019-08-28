"""
This file contains the complete code
to interact with the REACTOME Analysis
service's message queuing system.

This code makes use of the following environmental
parameters:
  * MAX_MESSAGE_TRIES
  * RABBIT_HOST
  * RABBIT_PORT
  * RABBIT_USER
  * RABBIT_PASSWORD
  * RABBIT_USER_SECRETS_FILE
  * RABBIT_PASSWORD_SECRETS_FILE
  * RABBIT_MAX_QUEUE_LENGTH
"""

import logging
import os
import signal

import pika
import pika.exceptions
import time

LOGGER = logging.getLogger(__name__)

# Names of queues to use for the different functions
ANALYSIS_QUEUE = "reactome_analysis_v0.1"
"""Queue used to send analysis requests to the worker"""
REPORT_QUEUE = "reactome_report_v0.1"
"""Queue used to send complete results to the report generating function"""
DATASET_QUEUE = "reactome_dataset_v0.1"
"""Queue used to retrieve external datasets"""


class ReactomeMQException(Exception):
    """
    Simple wrapper to identify exception created by the
    ReactomeMQ class.
    """
    pass


class ReactomeMQ:
    """
    Class used to manage all queuing systems in the Reactome Analysis System
    """

    """
    Currently, this class simply connects to the RabbitMQ
    instance to post to the analysis queue.
    """
    def __init__(self, queue_name=ANALYSIS_QUEUE):
        """Create a new ReactomeMQ instance

        :param queue_name: Name of the queue to use. This name should be taken from the queue names
                           in this package.
        """
        # use one queue per API version - in case the request format
        # changes, this needs to be changed as well
        self.queue_name = queue_name
        self.connection = None
        self._channel = None
        self._shutdown = False

        # stop analyses on TERM signals
        signal.signal(signal.SIGTERM, self._on_signal)
        signal.signal(signal.SIGINT, self._on_signal)

    def close(self):
        """
        Close the connection
        """
        if self.connection:
            self.connection.close()

    def post_analysis(self, analysis, method: str):
        """
        Post an analysis to the queue to be processed.
        :param analysis: The JSON-encoded analysis specification as a string.
        :param method: Name of the analysis method
        """
        # only allow a 3 second socket timeout for posting analysis requests
        try:
            channel = self._connect(3).channel()
        except Exception as e:
            raise ReactomeMQException("Failed to connect to queuing system: " + str(e))

        # declare the queue to make sure it exists
        channel.queue_declare(queue=self.queue_name, durable=True,
                              arguments={"x-overflow": "reject-publish",
                                         "x-max-length": int(os.getenv("RABBIT_MAX_QUEUE_LENGTH", 10))})

        # Turn on delivery confirmations
        channel.confirm_delivery()

        max_retries = int(os.getenv("MAX_MESSAGE_TRIES", 3))
        current_try = 1
        was_published = False

        # try to submit the job n times
        try:
            while current_try <= max_retries and not was_published:
                was_published = channel.basic_publish(exchange='',
                                                      routing_key=self.queue_name,
                                                      body=analysis,
                                                      properties=pika.BasicProperties(
                                                          delivery_mode=2  # make message persistent
                                                      ),
                                                      mandatory=True)  # require acknowledgement
                current_try += 1
                # wait for 1 sec to try again
                if not was_published:
                    LOGGER.debug("Failed to publish analysis message. Retrying...")
                    time.sleep(1)

            channel.close()
        except pika.exceptions.ConnectionClosed:
            raise ReactomeMQException("Message queue not accepting messages at the moment.")

        if not was_published:
            raise ReactomeMQException("Failed to publish analysis")

    def process_single_message(self, callback):
        """
        This function is mainly intended for testing purposes. It will connect, process a single message and return
        :param callback: Callback function to call
        """
        channel = self._connect().channel()
        method_frame, header_frame, body = channel.basic_get(self.queue_name)

        if body:
            callback(channel, method_frame, header_frame, body)

        channel.close()

    def sleep(self, duration: float) -> None:
        """
        "Sleep" for the given time. This should still send heartbeats
        while being a sleep.
        :param duration: The duration to sleep in seconds
        """
        self._connect().sleep(duration)

    def process_analysis(self, callback):
        """
        This function calls `callback` for every message it receives.

        The `callback` function must follow the format:
        ```
        def callback(ch, method, properties, body):
        ```
        and must acknowledge every message (once processed) with

        ```
        ch.basic_ack(delivery_tag=method.delivery_tag)
        ```

        This is a blocking function.

        :param callback: The function to call to process the message
        """
        while not self._shutdown:
            try:
                self._channel = self._connect().channel()

                self._channel.queue_declare(queue=self.queue_name, durable=True,
                                            arguments={"x-overflow": "reject-publish",
                                                       "x-max-length": int(os.getenv("RABBIT_MAX_QUEUE_LENGTH", 10))})

                # only allow one message to be posted at a time
                self._channel.basic_qos(prefetch_count=1)

                # register the callback
                self._channel.basic_consume(callback,
                                            queue=self.queue_name)

                # start listening
                LOGGER.info("Starting listening for analysis messages...")
                self._channel.start_consuming()
            except pika.connection.exceptions.ConnectionClosed:
                LOGGER.info("RabbitMQ connection closed")

                if not self._shutdown:
                    time.sleep(1)
                    LOGGER.info("Reconnecting...")
                    # simply re-connect
                    self.connection = None

                continue
            except pika.connection.exceptions.AMQPConnectionError:
                LOGGER.info("AMQPConnectionError received. Re-connecting...")
                # simply re-connect
                self.connection = None
                continue
            finally:
                self.close()

    def _on_signal(self, signum, frame):
        """
        Callback for SIGINT and SIGTERM signals
        :param signum:
        :param frame:
        """
        LOGGER.info("Signal receieved. Shutting down.")
        self.stop_analysis()

    def stop_analysis(self):
        """
        Stop the current analysis loop which was started
        using the `process_analysis` function.
        """
        LOGGER.info("Stopping analysis...")
        self._shutdown = True

        if self._channel:
            self._channel.stop_consuming()

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

    def _connect(self, timeout=pika.ConnectionParameters.DEFAULT_SOCKET_TIMEOUT):
        """
        Returns the active connection or create a new one.
        :param timeout: Timeout in (?seconds) for the socker connection.
        :return: A connection object
        """
        if self.connection is not None:
            return self.connection

        rabbit_host = os.getenv("RABBIT_HOST", "rabbit-mq")
        rabbit_port = os.getenv("RABBIT_PORT", 5672)

        # load user and password from the secrets file
        rabbit_user = "user"
        rabbit_password = "user"

        rabbit_file_user = ReactomeMQ._load_secret_from_file(os.getenv("RABBIT_USER_SECRETS_FILE", None))
        rabbit_file_pass = ReactomeMQ._load_secret_from_file(os.getenv("RABBIT_PASSWORD_SECRETS_FILE", None))

        if rabbit_file_user:
            rabbit_user = rabbit_file_user
        if rabbit_file_pass:
            rabbit_password = rabbit_file_pass

        # overwrite user and password with environmental values if set
        rabbit_env_user = os.getenv("RABBIT_USER", None)
        if rabbit_env_user:
            rabbit_user = rabbit_env_user

        rabbit_env_password = os.getenv("RABBIT_PASSWORD", None)
        if rabbit_env_password:
            rabbit_password = rabbit_env_password

        # connect to the service
        connection = pika.BlockingConnection(pika.ConnectionParameters(rabbit_host, port=rabbit_port,
                                                                       credentials=pika.PlainCredentials(
                                                                           username=rabbit_user,
                                                                           password=rabbit_password),
                                                                       socket_timeout=timeout,
                                                                       blocked_connection_timeout=timeout))

        # save the connection
        self.connection = connection

        return connection

    def get_is_shutdown(self):
        """
        Indicates whether the current analysis should be interrupted
        :return: boolean
        """
        return self._shutdown
