# system modules
import sys, os, shutil
from pathlib import Path
from contextlib import redirect_stdout
import multiprocessing

# text manipulation + configuration modules
import re
import linecache
import yaml
import logging

# numerical
import random
import numpy as np
import scipy.signal as sps
import pandas as pd

# plotting modules
import matplotlib.pyplot as plt
import scienceplots
import subprocess
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor

plt.style.use('science')
plt.rcParams["figure.figsize"] = (20,20)
plt.rcParams.update({
    'axes.labelsize': 14,      # Font size for axis labels
    'axes.titlesize': 16,      # Font size for axis titles
    'xtick.labelsize': 12,     # Font size for x-axis ticks
    'ytick.labelsize': 12,     # Font size for y-axis ticks
    #'font.family': 'serif',    # Font family (default is 'sans-serif')
    'font.size': 20,           # General font size for text in the plot
    'axes.labelsize': 20,
    'axes.titlesize': 24,
    'xtick.labelsize' : 16,
    'ytick.labelsize' : 16,
    'lines.linewidth':5,         # line width
    'lines.markersize' : 10,    # marker size
})

# constants
PI = 3.141592653589793
SBM_OUTPUT = "sbm_files"
SBM_DIR = os.path.dirname(os.path.realpath(__file__))

# functions
def truncate(val: float, stepsize: float) -> float:
    '''
    Truncates a value [val] based on the number of digits obtained from [stepsize].
    Example: val = 0.50000001; size = 0.01 => output = 0.5
    '''
    digit = 1.0/stepsize
    return round(val*digit) / digit

def get_openfoam_type(cfd_version: str) -> str:
    """
    Obtain the OpenFOAM type by checking the cfd_version string according to predefined regex patterns.
    E.g.    cfd_version = v2312-foss-2023a => "com"
            cfd_version = 11-foss-2023a => "org"

    Args:
        cfd_version: The OpenFOAM version.

    Returns:
        Returns Foundation (ORG) or ESI (COM)

    Raises:
        ValueError: Raises an exception if the cfd_version is unknown.
    """

    com_pattern = re.compile("^(v[0-9]+-foss-[0-9]+[a-b])$")
    com_match = re.match(com_pattern, cfd_version)
    if com_match:
        return "com"

    org_pattern = re.compile("^([0-9]+).*(-foss-).*([0-9]+[a-b])$")
    org_match = re.match(org_pattern, cfd_version)
    if org_match:
        return "org"

    raise ValueError("The cfd_version variable is unknown for OpenFOAM. \
                     Use the format [v2312-foss-2023a] for com or [11-foss-2023a] for org version")