"""
Copyright (C) 2023 Forschungszentrum Juelich GmbH, INM-6

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

"""
Writes out measures taken with internal instrumentation of the code.
MPI processes write to private logfiles. These files are scanned for the timer metrics.
Their mean is taken and writen into a single text file. This single text file will subsequently be read by JUBE.
"""

import glob
import os
import sys

import json
import yaml
import numpy as np
import pandas as pd


model_config_path = sys.argv[1]  # path to model config file
workpackage_path = sys.argv[2]   # path to workpackage files
log_path = sys.argv[3]           # path to logfiles


with open(model_config_path) as model_config_file:
    config = yaml.safe_load(model_config_file)

metrics = config["metrics"]
metrics_sum = config["metrics_sum"]

all_logfiles = glob.glob(os.path.join(log_path, '*logfile*'))

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

with open('timer_data.txt', "w") as out_file:
    for m in d:
        out_file.write(m + ' ' + str(d[m]) + '\n')
    for m in d_sum:
        out_file.write(m + ' ' + str(d_sum[m]) + '\n')

# CPU logging
cpu_info = [element.strip().replace(' ', '').replace('(', '').replace(')', '')
            for element in os.popen('lscpu').readlines()]

cpu_info_dict = {}
for element in cpu_info:
    key, value = element.split(':')
    cpu_info_dict[key] = value

with open(os.path.join(workpackage_path, 'cpu.json'), 'w') as f:
    json.dump(cpu_info_dict, f)

# Step timer and memory data
all_step_logfiles = glob.glob(os.path.join(log_path, '*logfile*_steps.dat'))

N = len(all_step_logfiles)
output_df = pd.read_csv(all_step_logfiles[0], sep=" ", header=None)
output_df.columns = ["time", "vmsize", "vmpeak", "vmrss"]
for logfile in all_step_logfiles[1:]:
    df = pd.read_csv(logfile, sep=" ", header=None)
    df.columns = ["time", "vmsize", "vmpeak", "vmrss"]
    output_df.add(df, fill_value=0)
output_df["time"] /= N
output_df["time_std"] = 0
for logfile in all_step_logfiles:
    df = pd.read_csv(logfile, sep=" ", header=None)
    output_df["time_std"] += (df["time"] - output_df["time"])**2
output_df["time_std"] /= N
output_df["time_std"] = np.sqrt(output_df["time_std"])

output_df.to_csv('timer_step_data.txt', sep='\t', encoding='utf-8', index=False)
