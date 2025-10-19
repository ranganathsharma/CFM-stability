"""
This file has the codes to calculate stability conditions for a given model at different traffic states and perturbation strengths.
The size of the platoon and duration of simulation are the hyperparameters of the system
"""

# Imports

import os, sys, yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, root_scalar
import logging
import multiprocessing as mp
from sa_idm import Stability_analysis

# Import the parameters

logging.basicConfig(filename= f'{__file__}.log', 
                    level=logging.WARNING, 
                    filemode = 'w')

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(parent_dir)

yaml_file_path = os.path.join(parent_dir, 'parameters.yaml')
with open(yaml_file_path, 'r') as file:
    model_param_dict = yaml.safe_load(file)

from model_equations import Model_equations
from solvers import solvers

print('The processing is complete')

model = 'idm'

dens = float(sys.argv[1]) # Input density from the command prompt

pert_kind = 'custom'

da_list = [0.1, 2, 5]
dtau_list = [1, 3, 5]

for da in da_list:
    for dtau in dtau_list:
        print(f'The process has begun for dens, da, dtau {dens, da, dtau} respectively')
        sa = Stability_analysis(model = model,
                                dens = dens,
                                pert_kind = pert_kind,
                                da = da,
                                dtau = dtau,
                                dt = 0.01)
        sa.looper()
        del sa
        print(f'Density, da, dtau = {dens, da, dtau} process completed')