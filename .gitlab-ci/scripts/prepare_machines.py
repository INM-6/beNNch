
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

# Remove temporary repositories again
shutil.rmtree(models_tmp_dir)

with open(variables_env_path, 'w') as variables_env:
    out_machines = list()
    for machine, machine_conf in config["user_config"].items():
        variables_env.write(f"CB_{machine.upper()}_PROJECT={machine_conf.get('project', '')}\n")
        variables_env.write(f"CB_{machine.upper()}_RUNNER_ACCOUNT={machine_conf.get('account', '')}\n")

    variables_env.write(f"CB_MACHINES={','.join(out_machines)}\n")
