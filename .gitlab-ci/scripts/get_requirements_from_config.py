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
Read all requirements from the configuration file and return them.
"""

import sys
import os

import yaml
yaml.add_representer(str, lambda dumper, data: dumper.represent_scalar('tag:yaml.org,2002:str', data, style='"'))


requirement_type = sys.argv[1]
config_file_path = sys.argv[2]
machine = sys.argv[3]

with open(config_file_path) as config_file:
    requirements = yaml.safe_load(config_file)["dependencies"][machine].get(requirement_type, [])

print(*requirements, sep=" ", end="")
