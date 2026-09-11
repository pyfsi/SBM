from utils import os, yaml

from src.sbm import SBM

if __name__=="__main__":
    config_path = os.path.join(os.getcwd(), "config.yaml")
    with open(config_path, "r") as conf_f:
        config = yaml.load(conf_f, Loader=yaml.SafeLoader)

    sbm = SBM(config)
    sbm.check_case()
    sbm.purge_boundary_data()
    sbm.initialize()
    sbm.run()
    sbm.finalize()
