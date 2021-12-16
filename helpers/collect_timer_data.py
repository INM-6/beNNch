"""
beNNch - Unified execution, collection, analysis and
comparison of neural network simulation benchmarks.
Copyright (C) 2021 Forschungszentrum Juelich GmbH, INM-6

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with
this program. If not, see <https://www.gnu.org/licenses/>.

SPDX-License-Identifier: GPL-3.0-or-later
"""

import glob
import os
import sys

import numpy as np

log_path = sys.argv[1]

"""
This function writes out measures taken with internal instrumentation of
the code. MPI processes write to private logfiles. These files are
scanned for the timer metrics. Their mean
is taken and writen into a single text file. This single text file can
later be read by eg JUBE.

Parameters
----------
STDOUT_PATH : string
    Place to store extracted timer data to
data_dir : string
    Directory of where simulation data is stored
label : string
    Unique identifier of a given simulation
"""

all_logfiles = glob.glob(
    os.path.join(
        log_path,
        '*logfile*'
    )
)

metrics = ['time_collocate_spike_data',
           'time_communicate_spike_data',
           'time_communicate_target_data',
           'time_deliver_spike_data',
           'time_gather_spike_data',
           'time_gather_target_data',
           'time_update',
           'time_communicate_prepare',
           'time_construction_connect',
           'time_construction_create',
           'time_simulate',
           'py_time_kernel_prepare',
           'py_time_network_local',
           'py_time_network_global',
           'py_time_simulate',
           'py_time_presimulate',
           'py_time_network_prepare',
           'py_time_create',
           'py_time_connect_area',
           'py_time_connect_cc',
           'py_time_connect']

metrics_sum = ['base_memory',
               'node_memory',
               'network_memory',
               'init_memory',
               'total_memory',
               'num_connections',
               'local_spike_counter']

d = {key: list() for key in metrics}
d_sum = {key: list() for key in metrics_sum}

for logfile in all_logfiles:
    with open(logfile, 'r') as fn:
        log = {}
        for line in fn:
            key, value, *_ = line.split(' ')
            if key in metrics + metrics_sum:
                log[key] = float(value)

        for m in d:
            try:
                d[m].append(log[m])
            except KeyError:
                pass
        for m in d_sum:
            try:
                d_sum[m].append(log[m])
            except KeyError:
                pass

for m in d:
    if d[m]:
        d[m] = np.mean(d[m])
    else:
        d[m] = np.nan

for m in d_sum:
    if d_sum[m]:
        d_sum[m] = np.sum(d_sum[m])
    else:
        d_sum[m] = np.nan

with open('timer_data.txt', "w") as outF:
    for m in d:
        outF.write(m + ' ' + str(d[m]) + '\n')
    for m in d_sum:
        outF.write(m + ' ' + str(d_sum[m]) + '\n')
