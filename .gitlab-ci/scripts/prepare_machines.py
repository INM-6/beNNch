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
import copy
import subprocess
import shutil

import yaml


config_file_path = sys.argv[1]
model_repos_dir = sys.argv[2]
machines = sys.argv[3]
variables_env_path = sys.argv[4]


with open(config_file_path) as config_file:
    config = yaml.safe_load(config_file)

combined_configs = dict()

models = os.environ.get("CB_MODELS", [])
if models:
    models = models.split(",")

models_tmp_dir = os.path.join(model_repos_dir, "tmp")
if not os.path.exists(models_tmp_dir):
    os.makedirs(models_tmp_dir)

model_tmp_dirs = dict()
for repo_alias, repo in config["repositories"].items():
    # Clone model repository into temporary location
    model_tmp_dirs[repo_alias] = os.path.join(models_tmp_dir, repo_alias)
    subprocess.run(["git", "clone", "--depth=1", "--single-branch", "--branch=" + repo["branch"], repo["url"], model_tmp_dirs[repo_alias]])

for model_name, benchmark_configs in config["benchmarks"].items():
    # Skip benchmarks defined in config which were not explicitly requested to be benchmarked
    if models and model_name not in models:
        continue
    
    # Copy relevant part of the repository to final destination
    model_repo_dir = os.path.join(model_repos_dir, model_name)
    shutil.copytree(os.path.join(model_tmp_dirs[benchmark_configs["repository"]], benchmark_configs["model_dir"]), model_repo_dir)
    # Copy .git for "model_git_commit_hash=$(git rev-parse HEAD)" in bennch helpers
    shutil.copytree(os.path.join(model_tmp_dirs[benchmark_configs["repository"]], ".git"), os.path.join(model_repo_dir, ".git"))

# Remove temporary repositories again
shutil.rmtree(models_tmp_dir)

with open(variables_env_path, 'w') as variables_env:
    out_machines = list()
    for machine, machine_conf in config["user_config"].items():
        out_machines.append(machine)
        variables_env.write(f"CB_{machine.upper()}_DATA_PROJECT={machine_conf.get('data-project', '')}\n")
        variables_env.write(f"CB_{machine.upper()}_COMPUTE_PROJECT={machine_conf.get('compute-project', '')}\n")
        variables_env.write(f"CB_{machine.upper()}_RUNNER_ACCOUNT={machine_conf.get('account', '')}\n")

    variables_env.write(f"CB_MACHINES={','.join(out_machines)}\n")
    variables_env.write(f"CB_DEPENDENCIES={config['dependencies']}")
