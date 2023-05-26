import sys
import yaml


config_file_path = sys.argv[1]
model_repos_dir = sys.argv[2]
models = sys.argv[3].split(",")


with open(config_file_path) as config_file:
    config = yaml.safe_load(config_file)

models_tmp_dir = os.path.join(model_repos_dir, "tmp")
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