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

import sys
import subprocess

import json
import numpy as np
import yaml
yaml.add_representer(str, lambda dumper, data: dumper.represent_scalar('tag:yaml.org,2002:str', data, style='"'))


result_file_path = sys.argv[1]
cpu_info = sys.argv[2]
config_file = sys.argv[3]
machine = sys.argv[4]
model_name = sys.argv[5]

with open(cpu_info) as f:
    for key, value in json.load(f).items():
        subprocess.run(['git', 'annex', 'metadata', result_file_path, '--set', f'{key}="{value.replace(" ", ";")}"'])

with open(config_file, 'r') as config_file:
    config = yaml.safe_load(config_file)

# TODO: Some parameters might not have concrete values, but instead either depend on other parameters or are based on python.
#  Instead of naively adding the value here, extract the actual values from jube somehow.
for k, v in config['benchmarks'][model_name]['machine_parameters'][machine].items():
    subprocess.run(['git', 'annex', 'metadata', result_file_path, '--set', f'{k}="{v}"'])

for v in config['benchmarks'][model_name]['model_parameters']:
    subprocess.run(['git', 'annex', 'metadata', result_file_path, '--set', f'{v["name"]}="{v["_"]}"'])
