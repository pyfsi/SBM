from utils import logging, os
from utils import contextmanager

class Logger():
    def __init__(self, config):
        self.logger = logging.getLogger(__name__)

        out_path = config.get("_output_path")
        log_path = os.path.join(out_path, "sbm.log")
        logging.basicConfig(filename=log_path, level=logging.INFO)

    def info(self, text):
        self.logger.info(text)

    def error(self, text):
        self.logger.error(text)

    @contextmanager
    def function_call(self, name):
        self.logger.info(f" ===== Start {name} ===== ")
        yield
        # TODO memory profiler
        self.logger.info(f" ===== End {name} ===== ")

    @contextmanager
    def subfunction_call(self, name):
        # TODO
        self.logger.info(f"")
        yield
        self.logger.info(f"")
