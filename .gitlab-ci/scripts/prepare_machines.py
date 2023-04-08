
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


for model_name, benchmark_configs in config["benchmarks"].items():
    # Skip benchmarks defined in config which were not explicitly requested to be benchmarked
    if models and model_name not in models:
        continue

    # Clone model repository into temporary location
    model_temp_dir = os.path.join(model_repos_dir, "tmp", model_name)
    subprocess.run(["git", "clone", "--depth=1", "--single-branch", "--branch=" + benchmark_configs["branch"], benchmark_configs["repository"], model_temp_dir])
    # Copy relevant part of the repository to final destination
    model_repo_dir = os.path.join(model_repos_dir, model_name)
    shutil.copytree(os.path.join(model_temp_dir, benchmark_configs["model_path"]), model_repo_dir)
    # Remove temporary repository again
    shutil.rmtree(model_temp_dir)

with open(variables_env_path, 'w') as variables_env:
    out_machines = list()
    for machine, machine_conf in config["user_config"].items():
        variables_env.write(f"CB_{machine.upper()}_PROJECT={machine_conf.get('project', '')}\n")
        variables_env.write(f"CB_{machine.upper()}_RUNNER_ACCOUNT={machine_conf.get('account', '')}\n")

    variables_env.write(f"CB_MACHINES={','.join(out_machines)}\n")
