import traceback
import unittest
import logging

class MyTestCase(unittest.TestCase):

    def get_func_name(self):
        stack = traceback.extract_stack()
        filename, codeline, funcName, text = stack[-3]
        return funcName

    def log_start(self):
        logging.debug('Start ' + self.get_func_name())