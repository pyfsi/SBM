import sys
import numpy as np
import os
import linecache
import random
import yaml
import shutil
from pathlib import Path
from contextlib import redirect_stdout
import multiprocessing
import logging
import re

from tqdm import tqdm
from time import sleep
import psutil

PI = 3.141592653589793
NPY_OUT_FOLDER = "sbm_files"

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