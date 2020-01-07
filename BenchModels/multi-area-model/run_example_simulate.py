"""
This script is used to run a simulation from the given command-line
arguments:
1. Scale of network
2. number of processes

It initializes the network class and then runs the simulate method of
the simulation class instance.

This script should be run after having run python run_example_createParams.py N_scaling num_processes sim_time K_scaling  as

mpirun/srun python run_example_simulate.py N_scaling num_processes nest_version
"""

import json
import nest
import os
import sys
import shutil

from config import base_path, data_path
from multiarea_model import MultiAreaModel, MultiAreaModel3

scale = float(sys.argv[1])
num_proc = int(sys.argv[2])
NEST_version = int(sys.argv[3])


# Load simulation and network labels
labels_fn = os.path.join(base_path, 'label_files/labels_{}_{}.json'.format(scale, num_proc))
with open(labels_fn, 'r') as f:
    labels = json.load(f)

label = labels['simulation_label']
network_label = labels['network_label']

# Load simulation parameters
fn = os.path.join(data_path, label, '_'.join(('custom_params', label)))
with open(fn, 'r') as f:
    custom_params = json.load(f)
# Copy custom param file for each MPI process
for i in range(custom_params['sim_params']['num_processes']):
    shutil.copy(fn, '_'.join((fn, str(i))))

fn = os.path.join(data_path,
                  label,
                  '_'.join(('custom_params',
                            label,
                            str(nest.Rank()))))
with open(fn, 'r') as f:
    custom_params = json.load(f)

os.remove(fn)

if NEST_version == 2:
    M = MultiAreaModel(network_label,
                       simulation=True,
                       sim_spec=custom_params['sim_params'])
elif NEST_version == 3:
    M = MultiAreaModel3(network_label,
                        simulation=True,
                        sim_spec=custom_params['sim_params'])

M.simulation.simulate()
