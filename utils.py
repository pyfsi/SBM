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

PI = 3.141592653589793
NPY_OUT_FOLDER = "sbm_files"

def truncate(val: float, stepsize: float) -> float:
    '''
    Truncates a value [val] based on the number of digits obtained from [stepsize].
    Example: val = 0.50000001; size = 0.01 => output = 0.5
    '''
    digit = 1.0/stepsize
    return round(val*digit) / digit