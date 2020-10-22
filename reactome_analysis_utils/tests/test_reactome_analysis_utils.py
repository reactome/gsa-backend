#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `reactome_analysis_utils` package."""


import unittest
import os
from reactome_analysis_utils.reactome_mq import ReactomeMQ

class ReactomeMqTest(unittest.TestCase):
    def setUp(self):
        os.environ["RABBIT_HOST"] = "127.0.0.1"
        os.environ["RABBIT_PORT"] = "30935"
        os.environ["RABBIT_USER"] = "test"
        os.environ["RABBIT_PASSWORD"] = "test"

    def test_connection(self):
        mq = ReactomeMQ()
        self.assertIsNotNone(mq)
        mq.close()

    def test_post_message(self):
        # disable posting using `rabbitmqctl set_vm_memory_high_watermark 0`
        mq = ReactomeMQ()
        mq.post_analysis("asda", "asda")
        mq.close()
