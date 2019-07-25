import logging

import connexion
import os
from flask_testing import TestCase

from reactome_analysis_api.encoder import JSONEncoder


class BaseTestCase(TestCase):

    def create_app(self):
        logging.getLogger('connexion.operation').setLevel('ERROR')
        app = connexion.App(__name__, specification_dir=os.path.join(os.path.dirname(__file__), "..",
                                                                     "reactome_analysis_api", "swagger"))
        app.app.json_encoder = JSONEncoder
        app.add_api('swagger.yaml')
        return app.app
