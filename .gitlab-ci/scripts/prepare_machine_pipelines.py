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


pipelines_dir_in = sys.argv[1]
pipelines_dir_out = sys.argv[2]
config_file_path = sys.argv[3]
model_repos_dir = sys.argv[4]
available_machines = sys.argv[5]
variables_env_path = sys.argv[6]


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


for machine, machine_conf in config["user_config"].items():
    with open(os.path.join(pipelines_dir_out, machine + ".yml"), 'w') as pipeline_out:
        with open(os.path.join(pipelines_dir_in, machine + ".yml")) as pipeline_in:
            pipeline_yml = yaml.safe_load(pipeline_in)

        if "variables" not in pipeline_yml:
            pipeline_yml["variables"] = dict()
        pipeline_yml["variables"]["CB_DATA_PROJECT"] = machine_conf.get('data-project', '')
        pipeline_yml["variables"]["CB_COMPUTE_PROJECT"] = machine_conf.get('compute-project', '')
        pipeline_yml["variables"]["CB_RUNNER_ACCOUNT"] = machine_conf.get('account', '')
        pipeline_yml["variables"]["CB_DEPENDENCIES"] = config['dependencies']
        yaml.dump(pipeline_yml, pipeline_out, default_flow_style=False)
        available_machines.remove(machine)


for machine in available_machines:
    print("Skipping pipeline generation for machine " + machine)
    with open(os.path.join(pipelines_dir_out, machine + ".yml"), 'w') as pipeline_out:
        pipeline_out.write(r"{stages: [fail], Fail: {stage: fail, rules: [{when: never}], script: ['']}}")
