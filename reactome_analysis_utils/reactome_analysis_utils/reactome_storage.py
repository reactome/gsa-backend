"""
Manages all functions required to connect to the
redis instance.

This script makes use of the following environmental variables (if present):

  * REDIS_HOST
  * REDIS_PORT
  * REDIS_PASSWORD_FILE: If present, this file will be used to load the password.
  * REDIS_PASSWORD: If present, this password will be use to connect to redis
                    and will overwrite a possible REDIS_PASSWORD_FILE content.
  * REDIS_DATABASE
"""

import logging
import os

import redis
import rediscluster

LOGGER = logging.getLogger(__name__)


class ReactomeStorageException(Exception):
    """
    Basic exception class to encapsulate issues
    related to the storage system.
    """
    pass


class ReactomeStorage:
    """
    Contains all functions to interact with the storage
    used by the REACTOME Analysis System
    """

    def __init__(self):
        # connect to redis
        try:
            self.r = self._get_redis()
        except Exception as e:
            raise ReactomeStorageException(e)

    def get_status(self, analysis_identifier: str, data_type: str = "analysis") -> str:
        """
        Returns the current status of the analysis.

        :param analysis_identifier: The analysis' identifier
        :param data_type: The data type to get the status for ["analysis", "report", "dataset"]
        :return: The status string (JSON-encoded) or None
        """
        try:
            if data_type == "report":
                status_key = self._get_report_status_key(analysis_identifier)
            elif data_type == "analysis":
                status_key = self._get_status_key(analysis_identifier)
            elif data_type == "dataset":
                status_key = self._get_request_data_status_key(analysis_identifier)
            else:
                raise ReactomeStorageException("Unknown type passed: " + data_type)

            LOGGER.debug("Getting status for {}".format(analysis_identifier))

            status = self.r.get(status_key)

            return status
        except Exception as e:
            raise ReactomeStorageException(e)

    def set_status(self, analysis_identifier: str, status: str, data_type: str = "analysis"):
        """
        Set the status of the analysis.

        :param analysis_identifier: The analysis' identifier
        :param data_type: The data type to get the status for ["analysis", "report", "dataset"]
        :param status: The status as JSON encoded string
        """
        try:
            if data_type == "report":
                status_key = self._get_report_status_key(analysis_identifier)
            elif data_type == "analysis":
                status_key = self._get_status_key(analysis_identifier)
            elif data_type == "dataset":
                status_key = self._get_request_data_status_key(analysis_identifier)
            else:
                raise ReactomeStorageException("Unknown type passed: " + data_type)

            LOGGER.debug("Setting status for {}: {}".format(analysis_identifier, status))
            self.r.set(status_key, status)
        except Exception as e:
            raise ReactomeStorageException(e)

    def get_result(self, analysis_identifier: str, data_type: str = "analysis") -> str:
        """
        Retrieve the result from storage as a JSON-encoded string.

        :param analysis_identifier: The analysis' identifier
        :param data_type: The data type to get the status for ["analysis", "report", "pdf_report]
        :return: The result as a string or None if it does not exist.
        """
        try:
            if data_type == "report":
                result_key = self._get_report_result_key(analysis_identifier)
            elif data_type == "pdf_report":
                result_key = self._get_pdf_report_result_key(analysis_identifier)
            elif data_type == "analysis":
                result_key = self._get_result_key(analysis_identifier)
            elif data_type == "r_script":
                result_key = self._get_r_script_result_key(analysis_identifier)
            else:
                raise ReactomeStorageException("Unknown type passed: " + data_type)

            LOGGER.debug("Getting result for {}".format(analysis_identifier))
            result = self.r.get(result_key)

            return result
        except Exception as e:
            raise ReactomeStorageException(e)

    def set_result(self, analysis_identifier: str, result: str, data_type: str = "analysis"):
        """
        Store the result for the specified analysis
        :param analysis_identifier: The analysis' identifier
        :param data_type: The data type to get the status for ["analysis", "report", "pdf_report", "r_script"]
        :param result: The result as a string
        """
        try:
            if data_type == "report":
                result_key = self._get_report_result_key(analysis_identifier)
            elif data_type == "pdf_report":
                result_key = self._get_pdf_report_result_key(analysis_identifier)
            elif data_type == "analysis":
                result_key = self._get_result_key(analysis_identifier)
            elif data_type == "r_script":
                result_key = self._get_r_script_result_key(analysis_identifier)
            else:
                raise ReactomeStorageException("Unknown type passed: " + data_type)

            LOGGER.debug("Setting result for {}".format(analysis_identifier))
            self.r.set(result_key, result)
        except Exception as e:
            raise ReactomeStorageException(e)

    def analysis_exists(self, analysis_identifier: str) -> bool:
        """
        Checks whether the specified analysis exists (has a status or result)
        :param analysis_identifier: The analysis' identifier
        :return: Boolean indicating whether the analysis exists
        """
        try:
            return self.r.exists(self._get_result_key(analysis_identifier)) or \
                   self.r.exists(self._get_status_key(analysis_identifier))
        except Exception as e:
            raise ReactomeStorageException(e)

    def set_request_data(self, token: str, data, expire: int = 3600):
        """
        Stores the passed request data under the specified token
        :param token: The token to store the data under
        :param data: The data to store
        :param expire: If not none, the key will be expired in `expire` seconds. Default = 60 Minutes = 3600 seconds.
        """
        try:
            request_key = self._get_request_data_key(token)

            self.r.set(request_key, data)

            if expire is not None and expire > 0:
                self.r.expire(request_key, expire)
        except Exception as e:
            raise ReactomeStorageException(e)

    def get_request_data(self, token: str) -> str:
        """
        Retrieve the stored request data for the given token
        :param token: The token under which the result was stored
        :return: The data
        """
        try:
            request_key = self._get_request_data_key(token)
            data = self.r.get(request_key)

            return data
        except Exception as e:
            raise ReactomeStorageException(e)

    def request_token_exists(self, token: str) -> bool:
        """
        Checks whether the request data token already exists.
        :param token: The token to check
        :return: Boolean indicating whether the token exists
        """
        try:
            return self.r.exists(self._get_request_data_key(token))
        except Exception as e:
            raise ReactomeStorageException(e)

    def set_request_data_summary(self, token: str, data, expire: int = 3600):
        """
        Stores the passed request data summary under the specified token
        :param token: The token to store the data under
        :param data: The data to store
        :param expire: If not none, the key will be expired in `expire` seconds. Default = 60 Minutes = 3600 seconds.
        """
        try:
            request_key = self._get_request_data_summary_key(token)

            self.r.set(request_key, data)

            if expire is not None and expire > 0:
                self.r.expire(request_key, expire)
        except Exception as e:
            raise ReactomeStorageException(e)

    def get_request_data_summary(self, token: str) -> str:
        """
        Retrieve the stored request data summary for the given token
        :param token: The token under which the result was stored
        :return: The data
        """
        try:
            request_key = self._get_request_data_summary_key(token)
            data = self.r.get(request_key)

            return data
        except Exception as e:
            raise ReactomeStorageException(e)

    def request_data_summary_exists(self, token: str) -> bool:
        """
        Checks whether the request data summary already exists.
        :param token: The token to check
        :return: Boolean indicating whether the token exists
        """
        try:
            return self.r.exists(self._get_request_data_summary_key(token))
        except Exception as e:
            raise ReactomeStorageException(e)

    def analysis_request_data_exists(self, token: str) -> bool:
        """
        Check whether the JSON-encoded analysis request object exists
        in storage.
        :param token: The token to check
        :return: Boolean indicating whether the token exists
        """
        try:
            return self.r.exists(self._get_analysis_request_key(token))
        except Exception as e:
            raise ReactomeStorageException(e)

    def get_analysis_request_data(self, token: str) -> str:
        """
        Retrieve the JSON-encoded analysis request object.
        :param token: The token identifying the analysis request.
        :return: The data
        """
        try:
            request_key = self._get_analysis_request_key(token)
            data = self.r.get(request_key)

            return data
        except Exception as e:
            raise ReactomeStorageException(e)

    def set_analysis_request_data(self, token: str, data:str, expire: int = 1800):
        """
        Store the JSON-encoded analysis request object.
        :param token: The token identifying the analysis request.
        :param data: The JSON-encoded string to store.
        :param expire: If not none, the key will be expired in `expire` seconds. Default = 30 Minutes = 1800 seconds.
        """
        try:
            request_key = self._get_analysis_request_key(token)
            
            self.r.set(name=request_key, value=data)

            if expire is not None and expire > 0:
                self.r.expire(request_key, expire)
        except Exception as e:
            raise ReactomeStorageException(e)

    @staticmethod
    def _get_redis():
        """
        Connects to the redis database as specified by the
        environmental variables.

        The redis API already uses an internal, thread-safe
        connection pool. Therefore, no additional connection
        pool has to be created.

        :return: A Redis connection pointing towards our redis instance.
        """

        # get the main config from environment variables or
        # use sensible default values
        redis_host = os.getenv("REDIS_HOST", "redis")
        redis_port = int(os.getenv("REDIS_PORT", 6379))
        redis_database = int(os.getenv("REDIS_DATABASE", 0))
        use_redis_cluster = os.getenv("USE_REDIS_CLUSTER", "False") == "True"

        # load the password from file if set
        redis_password_file = os.getenv("REDIS_PASSWORD_FILE", None)
        redis_password = None

        if redis_password_file and os.path.isfile(redis_password_file):
            with open(redis_password_file, "r") as reader:
                for line in reader:
                    line = line.strip()
                    if len(line) > 0:
                        redis_password = line

        # environment variable overwrites file
        redis_env_password = os.getenv("REDIS_PASSWORD", None)

        if redis_env_password:
            redis_password = redis_env_password

        LOGGER.debug("Connecting to Redis at " + redis_host + ":" + str(redis_port))

        # TODO: select redis cluster or the default redis depending on some (new) config variable
        startup_nodes = [{"host": redis_host, "port": redis_port}]

        if use_redis_cluster:
            redis_connection = rediscluster.RedisCluster(startup_nodes=startup_nodes, 
                                                         password=redis_password,
                                                         skip_full_coverage_check=True,
                                                         retry_on_timeout=False, 
                                                         socket_keepalive=False, 
                                                         socket_timeout=3,
                                                         socket_connect_timeout=3)
        else:
            redis_connection = redis.Redis(host=redis_host, port=redis_port, db=redis_database, password=redis_password,
                                           retry_on_timeout=False, socket_keepalive=False, socket_timeout=3,
                                           socket_connect_timeout=3)

            return redis_connection

    @staticmethod
    def _get_result_key(analysis_id: str) -> str:
        """
        Creates the redis result key for the specified analysis id.

        :param analysis_id: The analysis id.
        :return: The matching redis result key
        """
        return "analysis:{}:result".format(analysis_id)

    @staticmethod
    def _get_report_result_key(token: str) -> str:
        """
        Creates the redis key for the specified report data

        :param token: The token identifying the report
        :return: The matching redis key
        """
        return "report:{}:result".format(token)

    @staticmethod
    def _get_pdf_report_result_key(token: str) -> str:
        """
        Creates the redis key for the specified PDF report data

        :param token: The token identifying the report
        :return: The matching redis key
        """
        return "pdf_report:{}:result".format(token)

    @staticmethod
    def _get_r_script_result_key(token: str) -> str:
        """
        Creates the redis key for the specified R report script

        :param token: The token identifying the report
        :return: The matching redis key
        """
        return "r_script:{}:result".format(token)

    @staticmethod
    def _get_status_key(analysis_id: str) -> str:
        """
        Creates the redis status key for the specified analysis id.

        :param analysis_id: The analysis id.
        :return: The matching redis status key
        """
        return "analysis:{}:status".format(analysis_id)

    @staticmethod
    def _get_report_status_key(report_id: str) -> str:
        """
        Creates the redis status key for the specified report.

        :param report_id: The report id
        :return: The matching redis report id
        """
        return "report:{}:status".format(report_id)

    @staticmethod
    def _get_dataset_status_key(dataset_id: str) -> str:
        """
        Creates the redis status key for the specified dataset status.

        :param dataset_id: The dataset id
        :return: The matching redis report id
        """
        return "dataset:{}:status".format(dataset_id)

    @staticmethod
    def _get_request_data_key(token: str) -> str:
        """
        Creates the redis key for the specified request data

        :param token: The token identifying the request data
        :return: The matching redis key
        """
        return "request_data:{}".format(token)

    @staticmethod
    def _get_request_data_status_key(token: str) -> str:
        """
        Creates the redis key for the specified request data

        :param token: The token identifying the request data
        :return: The matching redis key
        """
        return "request_data:{}:status".format(token)

    @staticmethod
    def _get_request_data_summary_key(token: str) -> str:
        """
        Creates the redis key for the specified request data summary

        :param token: The token identifying the request data
        :return: The matching redis key
        """
        return "request_data:{}:summary".format(token)

    @staticmethod
    def _get_analysis_request_key(token: str) -> str:
        """
        Creates the redis key for the AanalysisRequest object (stored as JSON)

        :param token: The token identifying the analysis request
        :return: The matching redis key
        """
        return "analysis_request:{}:data".format(token)
