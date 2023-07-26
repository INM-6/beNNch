import sys
import os
import subprocess

import yaml


tmp_dir = sys.argv[1]
config_out_path = sys.argv[2]
config_repo = sys.argv[3]
config_branch = sys.argv[4]
config_path = sys.argv[5]


def patch_element(config, key, value):
    if key not in config:
        config[key] = value
    else:
        if type(config[key]) != type(value):
            raise ValueError()
        if type(config[key]) == list:
            config[key].sort()
            value.sort()
            for 
        if type(config[key]) == dict:
            for k, v in value.items():
                patch_element(config[key], k, v)
        else:
            config[key] = value


def extend_config(config, level):
    base_repo_alias = config["extends"]["repo"]
    base_repo_url = config["repositories"][base_repo_alias]["url"]
    base_repo_branch = config["repositories"][base_repo_alias]["branch"]
    base_repo_config_path = config["extends"]["config_file"]

    repo_path = os.path.join(tmp_dir, str(level))
    subprocess.run(["git", "clone", "--depth=1", "--single-branch", "--branch=" + base_repo_branch, base_repo_url, repo_path])
    with open(os.path.join(repo_path, base_repo_config_path), 'r') as config_file:
        base_config = yaml.safe_load(config_file)
    # Patch dependencies
    base_config["dependencies"] = os.path.join(repo_path, base_config["dependencies"])

    if "extends" in base_config:
        base_config = extend_config(base_config, level+1)

    for key, value in config.items():
        patch_element(base_config, key, value)
    return base_config

# Recursively get extended config file until reaching root. Starting with the root config, iteratively patch the configuration.
config = extend_config({"extends": {"repo": "base", "config_file": config_path}, "repositories": {"base": {"url": config_repo, "branch": config_branch}}}, 0)

# Load dependencies
with open(config["dependencies"], 'r') as dependencies_file:
    base_config["dependencies"] = yaml.safe_load(dependencies_file)

# TODO: Test this
# After patching the configuration using various base configurations, patch the final config using environment variables provided by the user when triggering the pipeline.
# for env_var_key, env_var_value in os.environ.items():
#     if env_var_key.startswith("CB_"):
#         keys = env_var_key[3:].split("__")
#         try:
#             reduce(operator.getitem, keys[:-1], config)[keys[-1]] = env_var_value
#         except:
#             print(f"The configuration path {'.'.join(keys[:-1])} is not present in original config. The variable {env_var_key} with value {env_var_value} can therefore not be set!", file=sys.stderr)


# Write patched config file
with open(config_out_path, 'w') as config_outfile:
    yaml.dump(config, config_outfile, default_flow_style=False)
