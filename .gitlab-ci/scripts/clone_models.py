import sys
import os
import subprocess
import copy
import shutil
from ruamel.yaml import YAML
yaml = YAML(typ='safe')
yaml.preserve_quotes = True
yaml.default_style = '"'


config_file_path = sys.argv[1]
model_repos_dir = sys.argv[2]
models = sys.argv[3]
if models:
    models = models.split(",")


with open(config_file_path) as config_file:
    config = yaml.load(config_file)

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
        repo_path = model_tmp_dirs[repo_alias]
        os.makedirs(repo_path, exist_ok=False)
        subprocess.run(["git", "init", "-q"],  cwd=repo_path, check=True, stdout=subprocess.DEVNULL)
        subprocess.run(["git", "fetch", "-q", "--depth=1", repo["url"], repo["branch"]], cwd=repo_path, check=True, stdout=subprocess.DEVNULL)
        subprocess.run(["git", "checkout", "-qf", "FETCH_HEAD"], cwd=repo_path, check=True, stdout=subprocess.DEVNULL)

    # Copy relevant part of the repository to final destination
    model_repo_dir = os.path.join(model_repos_dir, model_name)
    shutil.copytree(os.path.join(model_tmp_dirs[repo_alias], benchmark_configs["model_dir"]), model_repo_dir)
    # Copy .git for "model_git_commit_hash=$(git rev-parse HEAD)" in bennch helpers
    shutil.copytree(os.path.join(model_tmp_dirs[repo_alias], ".git"), os.path.join(model_repo_dir, ".git"), dirs_exist_ok=True)

# Remove temporary repositories again
shutil.rmtree(models_tmp_dir)
