from utils import logging, os
from utils import contextmanager
from utils import tracemalloc, time
from utils import is_inside

ONE_KIBIBYTE = 1 << 10
ONE_MEBIBYTE = 1 << 20
ONE_GIBIBYTE = 1 << 30
mem_unit_size = ONE_MEBIBYTE

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
        with self.memory_profiler(name):
            with self.time_profiler(name):
                yield
        self.logger.info(f" ===== End {name} ===== ")

    @contextmanager
    def subfunction_call(self, name):
        self.logger.info(f"")
        yield
        self.logger.info(f"")

    @contextmanager
    def time_profiler(self, name: str):
        start_time = time.time()
        yield
        run_time = time.time() - start_time
        self.logger.info(f"{name} runtime = {run_time:,.3f} seconds")

    @contextmanager
    def memory_profiler(self, name: str):
        tracemalloc.start()
        yield
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        # choose memory unit for printing
        if is_inside(peak, min=ONE_GIBIBYTE, max=100*ONE_GIBIBYTE):
            mem_unit = "GiB"
            mem_unit_size = ONE_GIBIBYTE
        elif is_inside(peak, min=ONE_MEBIBYTE, max=ONE_GIBIBYTE):
            mem_unit = "MiB"
            mem_unit_size = ONE_MEBIBYTE
        elif is_inside(peak, min=0, max=ONE_MEBIBYTE):
            mem_unit = "KiB"
            mem_unit_size = ONE_KIBIBYTE
        else:
            mem_unit = "GiB"
            mem_unit_size = ONE_GIBIBYTE

        self.logger.info(f"{name} memory allocation [{mem_unit}]: Peak = {peak / mem_unit_size:,.3f}; Final = {current / mem_unit_size:,.3f}")