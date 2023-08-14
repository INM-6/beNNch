import sys
import os
import subprocess
from functools import reduce
import operator

import yaml
yaml.add_representer(str, lambda dumper, data: dumper.represent_scalar('tag:yaml.org,2002:str', data, style='"'))

COLOR_WARNING = '\033[93m'
COLOR_FAIL = '\033[91m'
COLOR_END = '\033[0m'

tmp_dir = sys.argv[1]
config_out_path = sys.argv[2]
config_repo = sys.argv[3]
config_branch = sys.argv[4]
config_path = sys.argv[5]


def patch_element(config, key, value):
    if key not in config:
        config[key] = value
        print(f"patching key {key}", file=sys.stderr)
    else:
        print(f"patching key {key} with value of type {type(config[key])}", file=sys.stderr)
        if type(config[key]) != type(value):
            raise ValueError()
        if type(config[key]) == list:
            for v in value:
                if type(v) == dict and "name" in v:
                    print(v, file=sys.stderr)
                    for val in config[key]:
                        print("val:", val, file=sys.stderr)
                        if val.get("name", "") == v["name"]:
                            if "type" in v:
                                val["type"] = v["type"]
                            val["_"] = v["_"]
                    else:
                        config[key].append(v)
                elif v not in config[key]:
                    config[key].append(v)
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
    subprocess.run(["git", "clone", "-q", "--depth=1", "--single-branch", "--branch=" + base_repo_branch, base_repo_url, repo_path], stdout=subprocess.DEVNULL)
    with open(os.path.join(repo_path, base_repo_config_path), 'r') as config_file:
        base_config = yaml.safe_load(config_file)
    # Patch dependencies
    if "dependencies" in base_config:
        base_config["dependencies"] = os.path.join(repo_path, base_config["dependencies"])

    if "extends" in base_config:
        base_config = extend_config(base_config, level+1)

    for key, value in config.items():
        patch_element(base_config, key, value)
    print(base_config)
    return base_config

# Recursively get extended config file until reaching root. Starting with the root config, iteratively patch the configuration.
config = extend_config({"extends": {"repo": "base", "config_file": config_path}, "repositories": {"base": {"url": config_repo, "branch": config_branch}}}, 0)

# Load dependencies
with open(config["dependencies"], 'r') as dependencies_file:
    config["dependencies"] = yaml.safe_load(dependencies_file)

# After patching the configuration using various base configurations, patch the final config using environment variables provided by the user when triggering the pipeline.
for env_var_key, env_var_value in os.environ.items():
    if env_var_key.startswith("CB_") and env_var_key not in ("CB_MODELS", "CB_MACHINES"):
        keys = env_var_key[3:].split("__")
        try:
            reduce(operator.getitem, keys[:-1], config)[keys[-1]] = env_var_value
        except Exception as e:
            print(f"{COLOR_WARNING}The configuration path {'.'.join(keys[:-1])} is not present in original config. The variable {env_var_key} with value {env_var_value} can therefore not be set!{COLOR_END}", file=sys.stderr)

success = False
models = os.environ.get("CB_MODELS", "")
if models:
    models = models.split(",")

if os.environ["CB_MACHINES"]:
    machines = list(set(os.environ["CB_MACHINES"].split(",")) & set(config["user_config"].keys()))
else:
    machines = list(config["user_config"].keys())

machine_models = {machine: 0 for machine in machines}
to_remove = list()
for model_name, benchmark_config in config["benchmarks"].items():
    if model_name not in models:
        to_remove.append(model_name)
        continue

    if models and model_name not in models:
        print(f"{COLOR_WARNING}Skipping config generation for model {model_name} as requested by the user!{COLOR_END}", file=sys.stderr)
        continue
    for machine in list(benchmark_config["machine_parameters"].keys()):
        if machine not in machines:
            del benchmark_config["machine_parameters"][machine]
        else:
            machine_models[machine] += 1
for model in to_remove:
    del config["benchmarks"][model]

for machine, count in machine_models.items():
    if count == 0:
        machines.remove(machine)

# Write patched config file
with open(config_out_path, 'w') as config_outfile:
    yaml.dump(config, config_outfile, default_flow_style=False)

print(",".join([machine for machine, count in machine_models.items() if count != 0]))
