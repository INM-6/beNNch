import sys
import os
import subprocess
import copy
import shutil
import yaml
yaml.add_representer(str, lambda dumper, data: dumper.represent_scalar("tag:yaml.org,2002:str", data, style='"'))


config_file_path = sys.argv[1]
model_repos_dir = sys.argv[2]
models = sys.argv[3]
if models:
    models = models.split(",")


with open(config_file_path) as config_file:
    config = yaml.safe_load(config_file)

models_tmp_dir = os.path.join(model_repos_dir, "tmp")
os.makedirs(models_tmp_dir)

model_tmp_dirs = dict()

for model_name, benchmark_configs in config["benchmarks"].items():
    # Skip benchmarks defined in config which were not explicitly requested to be benchmarked
    if models and model_name not in models:
        continue
    
    repo_alias = benchmark_configs["repository"]
    if repo_alias not in model_tmp_dirs:
        model_tmp_dirs[repo_alias] = os.path.join(models_tmp_dir, repo_alias)
        repo = config["repositories"][repo_alias]
        subprocess.run(["git", "clone", "--depth=1", "--single-branch", "--branch=" + repo["branch"], repo["url"], model_tmp_dirs[repo_alias]])

    # Copy relevant part of the repository to final destination
    model_repo_dir = os.path.join(model_repos_dir, model_name)
    shutil.copytree(os.path.join(model_tmp_dirs[repo_alias], benchmark_configs["model_dir"]), model_repo_dir)
    # Copy .git for "model_git_commit_hash=$(git rev-parse HEAD)" in bennch helpers
    shutil.copytree(os.path.join(model_tmp_dirs[repo_alias], ".git"), os.path.join(model_repo_dir, ".git"))

# Remove temporary repositories again
shutil.rmtree(models_tmp_dir)

# Clone the simulator
env = os.environ.copy()
env["CB_SIMULATOR_REPO"] = os.path.expandvars(config["software_parameters"]["simulator_repo_url"])
#os.environ["CB_SIMULATOR_REPO"] = config["software_parameters"]["simulator_repo_url"]
simulator_branch = os.path.expandvars(config["software_parameters"]["simulator_repo_branch"])
subprocess.run(["build", "--get", "--silent", "simulator", simulator_branch], env=env, stdout=subprocess.DEVNULL)
