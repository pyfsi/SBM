# Synthetic Bubble Model

This repository contains the complete model used to define a transient inlet boundary condition to be used in a subsequent Volume-Of-Fluid (VOF) simulation. This model was coined the "Synthetic Bubble Model" (SBM).

The model constructs a virtual pre-domain in front of the actual domain's inlet and subsequently fills the pre-domain - initially filled with liquid - with gas bubbles of arbitrary shape and of random size, at a randomly chosen position in time and space. It is initially constructed for use in combination with OpenFOAM/4.1, but has been tested in OpenFOAM/6 as well and has been adapted for Ansys Fluent. Moreover, the main modelling script is independent of flow solver to be used in the subsequent CFD simulation.

The complete user manual "inletmodelling_manual" is also given in this repository; this a technical manual, without explanation of the concept behind the model or the underlying algorithm. The latter are described in the Ph.D. dissertation of Laurent De Moerloose, published in 2020 at Ghent University. You can find the full text [here](https://lib.ugent.be/catalog/rug01:002978914). The user manual is also in the appendix of this dissertation.

# How To Run
The configuration settings for the SBM is stored in `config.yaml`. The absolute path to the CFD case can be defined in the `case_path` variable inside `config.yaml`. An example configuration file can be found in this repository.

The following lines of code can be run in the terminal in order to run the code:
```
export PYTHON_VERSION=Anaconda3-python/2023.09-0
module load $PYTHON_VERSION
PYTHONDONTWRITEBYTECODE=1 python main.py -p no:cacheprovider
```
The specified Python version from an Anaconda distribution is loaded using the `module load` command. Then, the script `main.py` can be run by typing `python main.py`. This will create `__pycache__` folders containing bytecode-compiled versions of the program, which makes the program run a little bit faster. The generation of these `__pycache__` folders can be deactivated by running `PYTHONDONTWRITEBYTECODE=1 python main.py -p no:cacheprovider`. Additionally, this program can be executed by running the `master.sh` script using the `bash` command.

