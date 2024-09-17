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

user_machine = os.environ["CB_MACHINES"]

with open(config_file_path) as config_file:
    config = yaml.safe_load(config_file)

machine_conf = config["user_config"][user_machine]
print(machine_conf)
with open(os.path.join(pipelines_dir_in)) as pipeline_in:
          pipeline_yml = yaml.safe_load(pipeline_in)

pipeline_yml["variables"]["CB_DATA_PROJECT"] = machine_conf.get('data-project', '')
pipeline_yml["variables"]["CB_COMPUTE_PROJECT"] = machine_conf.get('compute-project', '')
pipeline_yml["variables"]["CB_RUNNER_ACCOUNT"] = machine_conf.get('account', '')
pipeline_yml["variables"]["CB_BUILD_TIME_LIMIT"] = machine_conf.get('build-time-limit', '')
pipeline_yml["variables"].update(config["dependencies"][user_machine]["variables"])

# get the desired simulator version
pipeline_yml["variables"]["CB_SIMULATOR_REPO"] = config["software_parameters"]["simulator_repo_url"]
pipeline_yml["variables"]["CB_SIMULATOR_BRANCH"] = config["software_parameters"]["simulator_repo_branch"]

# write package dependencies to the machine config
pipeline_yml["variables"]["PIP_LIST"] = config["dependencies"][user_machine]["variables"]["PIP_LIST"]
pipeline_yml["variables"]["MODULE_LIST"] = config["dependencies"][user_machine]["variables"]["MODULE_LIST"]

benchmarks = list()
benchmark_id = 0
for model_name, model_config in config["benchmarks"].items():
    if user_machine in model_config["machine_parameters"]:
        benchmarks.append({"JUBE_BENCHMARK_FILE": model_name,
               		   "JUBE_BENCHMARK_ID": str(benchmark_id).zfill(6)})
        benchmark_id += 1

pipeline_yml[".parallel-benchmark"] = {"parallel": {"matrix": benchmarks}}
print(os.path.join(pipelines_dir_out, user_machine + ".yml"))
with open(os.path.join(pipelines_dir_out), 'w') as pipeline_out:
          yaml.dump(pipeline_yml, pipeline_out, default_flow_style=False)

