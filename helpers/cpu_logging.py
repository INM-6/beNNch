"""
NEST Benchmarking Framework - Unified execution, collection, analysis and
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

import os
import sys
import json

save_path = sys.argv[1]

cpu_info = [element.strip().replace(' ', '').replace('(', '').replace(')', '')
            for element in os.popen('lscpu').readlines()]

cpu_info_dict = {}
for element in cpu_info:
    key, value = element.split(':')
    cpu_info_dict[key] = value

with open(os.path.join(save_path, 'cpu.json'), 'w') as f:
    json.dump(cpu_info_dict, f)
