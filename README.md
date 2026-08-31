# Synthetic Bubble Model

This repository contains the complete model used to define a transient inlet boundary condition to be used in a subsequent Volume-Of-Fluid (VOF) simulation. This model was coined the "Synthetic Bubble Model" (SBM).

The model constructs a virtual pre-domain in front of the actual domain's inlet and subsequently fills the pre-domain - initially filled with liquid - with gas bubbles of arbitrary shape and size at a randomly chosen position in time and space. It is initially constructed for use in combination with OpenFOAM/4.1, but has been tested in OpenFOAM/6 as well and has been adapted for Ansys Fluent. Moreover, the main modelling script is independent of the CFD flow solver.

The complete user manual "inletmodelling_manual" is also given in this repository; this a technical manual, without explanation of the concept behind the model or the underlying algorithm. The latter are described in the Ph.D. dissertation of Laurent De Moerloose, published in 2020 at Ghent University. You can find the full text [here](https://lib.ugent.be/catalog/rug01:002978914). The user manual is also in the appendix of this dissertation.

# How To Run
In Linux, the bash script `setup.sh` is used to setup the necessary environment variables. This bash script creates an alias to the SBM directory and adds it to `PYTHONPATH` as well as loads the specified Python version. Alternatively, the user can utilize their own Python installation and add the following line containing the SBM directory path to their `.bashrc` file:
```
export SBM=/path/to/sbm/directory
export PYTHONPATH=/path/to/sbm/directory:$PYTHONPATH
```

The configuration settings for the SBM is stored in the `config.yaml` file located inside the case directory, where the CFD simulation files are found. An example configuration file can be found in the `examples/channel_bubble` directory.

The following lines of code can be executed in the terminal in order to run the code:
```
python $SBM/main.py
```
Executing the script in this manner will create `__pycache__` files which improves execution time during subsequent runs. It is also possible to prevent Python from creating `__pycache__` files by running the script in the following way:
```
PYTHONDONTWRITEBYTECODE=1 python $SBM/main.py -p no:cacheprovider
```
