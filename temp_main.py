from utils import os, yaml, shutil
from utils import get_openfoam_type, memory_profiler, modulo
from utils import SBM_OUTPUT
import cProfile

from src_oop.sbm import SBM

if __name__=="__main__":
    config_path = os.path.join(os.getcwd(), "config.yaml")
    with open(config_path, "r") as conf_f:
        config = yaml.load(conf_f, Loader=yaml.SafeLoader)

    # SBM object
    sbm = SBM(config)
    sbm.check_case()
    sbm.purge_previous()
    sbm.initialize()
    sbm.run()
    sbm.finalize()
