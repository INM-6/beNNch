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
Scan the configuration file for all repositories and clone them.
Collect all machines for which benchmarks should be run and save them to an environment variable.
"""

import sys
import os
import traceback

import yaml
yaml.add_representer(str, lambda dumper, data: dumper.represent_scalar('tag:yaml.org,2002:str', data, style='"'))


pipelines_dir_in = sys.argv[1]
pipelines_dir_out = sys.argv[2]
config_file_path = sys.argv[3]
available_machines = sys.argv[4].split(",")

manual = True
if os.environ["CB_MACHINES"]:
    available_machines = os.environ["CB_MACHINES"].split(",")
    manual = False


with open(config_file_path) as config_file:
    config = yaml.safe_load(config_file)

for machine, machine_conf in config["user_config"].items():
    try:
        with open(os.path.join(pipelines_dir_in, machine + ".yml")) as pipeline_in:
            pipeline_yml = yaml.safe_load(pipeline_in)

        if "variables" not in pipeline_yml:
            pipeline_yml["variables"] = dict()
        if ".environment-init" not in pipeline_yml:
            pipeline_yml[".environment-init"] = list()
        if ".environment-load" not in pipeline_yml:
            pipeline_yml[".environment-load"] = list()
        if "variables" not in pipeline_yml:
            pipeline_yml["variables"] = list()
        if ".preparation" not in pipeline_yml:
            pipeline_yml[".preparation"] = dict()
        pipeline_yml["variables"]["CB_DATA_PROJECT"] = machine_conf.get('data-project', '')
        pipeline_yml["variables"]["CB_COMPUTE_PROJECT"] = machine_conf.get('compute-project', '')
        pipeline_yml["variables"]["CB_RUNNER_ACCOUNT"] = machine_conf.get('account', '')
        pipeline_yml["variables"]["CB_BUILD_TIME_LIMIT"] = machine_conf.get('build-time-limit', '')
        pipeline_yml["variables"].update(config["dependencies"][machine]["variables"])

        pipeline_yml[".environment-init"].extend(config["dependencies"][machine]["init"])
        pipeline_yml[".environment-load"].extend(config["dependencies"][machine]["load"])

        if manual:
            pipeline_yml[".preparation"]["when"] = "manual"
        else:
            pipeline_yml[".preparation"]["when"] = "always"
    except Exception as e:
        traceback.print_exc()
        continue

    with open(os.path.join(pipelines_dir_out, machine + ".yml"), 'w') as pipeline_out:
        yaml.dump(pipeline_yml, pipeline_out, default_flow_style=False)
        available_machines.remove(machine)

# All machines for which no config exists are ignored
for machine in available_machines:
    print("Skipping pipeline generation for machine " + machine)
    with open(os.path.join(pipelines_dir_out, machine + ".yml"), 'w') as pipeline_out:
        pipeline_out.write(r"{stages: [fail], Fail: {stage: fail, rules: [{when: never}], script: ['']}}")
