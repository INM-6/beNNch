import sys
import os
import shutil
import subprocess
from functools import reduce
import operator

from ruamel.yaml import YAML
yaml = YAML(typ='safe')
yaml.default_flow_style = False
yaml.preserve_quotes = True
yaml.default_style = '"'

COLOR_WARNING = '\033[93m'
COLOR_FAIL = '\033[91m'
COLOR_END = '\033[0m'

tmp_dir = sys.argv[1]
config_out_path = sys.argv[2]
config_repo = sys.argv[3]
config_branch = sys.argv[4]
config_path = sys.argv[5]


def lower_keys(x):
    if isinstance(x, list) or isinstance(x, tuple):
        return [lower_keys(v) for v in x]
    elif isinstance(x, dict):
        return {k.lower(): lower_keys(v) for k, v in x.items()}
    else:
        return x


def patch_element(config, key, value):
    if key not in config:
        config[key] = value
    else:
        if type(config[key]) != type(value):
            raise ValueError()
        if type(value) == list:
            if len(value) > 0 and type(value[0]) == dict:
                for v in value:
                    assert type(v) == dict, "Configuration files do not allow lists with mixed child types (e.g. scalar and dict)"
                    if "name" in v:
                        for i, val in enumerate(config[key]):
                            if val.get("name", "") == v["name"]:
                                config[key][i] = v
                                break
                        else:
                            config[key].append(v)
            else:
                config[key] = value
        elif type(value) == dict:
            for k, v in value.items():
                patch_element(config[key], k, v)
        else:
            config[key] = value


def extend_config(config, level):
    base_repo_alias = config["extends"]["repo"]
    base_repo_url = config["repositories"][base_repo_alias]["url"].replace("${CI_JOB_TOKEN}", os.environ['CI_JOB_TOKEN'])
    base_repo_branch = config["repositories"][base_repo_alias]["branch"]
    base_repo_config_path = config["extends"]["config_file"]

    repo_path = os.path.join(tmp_dir, str(level))
    os.makedirs(repo_path, exist_ok=True)
    subprocess.run(["git", "init", "-q"],  cwd=repo_path, check=True, stdout=subprocess.DEVNULL)
    subprocess.run(["git", "fetch", "-q", "--depth=1", base_repo_url, base_repo_branch], cwd=repo_path, check=True, stdout=subprocess.DEVNULL)
    subprocess.run(["git", "checkout", "-qf", "FETCH_HEAD"], cwd=repo_path, check=True, stdout=subprocess.DEVNULL)

    with open(os.path.join(repo_path, base_repo_config_path), 'r') as config_file:
        base_config = yaml.load(config_file)
    # Patch dependencies
    if "dependencies" in base_config:
        base_config["dependencies"] = os.path.join(repo_path, base_config["dependencies"])

    if "extends" in base_config:
        base_config = extend_config(base_config, level+1)
    for key, value in config.items():
        patch_element(base_config, key, value)
    return base_config

# Recursively get extended config file until reaching root. Starting with the root config, iteratively patch the configuration.
config = lower_keys(extend_config({"extends": {"repo": "base", "config_file": config_path}, "repositories": {"base": {"url": config_repo, "branch": config_branch}}}, 0))

# Load dependencies
with open(config["dependencies"], 'r') as dependencies_file:
    config["dependencies"] = yaml.load(dependencies_file)

# After patching the configuration using various base configurations, patch the final config using environment variables provided by the user when triggering the pipeline.
for env_var_key, env_var_value in os.environ.items():
    if env_var_key.startswith("CB_") and env_var_key.lower() not in ("cb_models", "cb_machines"):
        keys = [k.lower() for k in env_var_key[3:].split("__")]
        try:
            reduce(operator.getitem, keys[:-1], config)[keys[-1]] = env_var_value
        except Exception as e:
            print(f"{COLOR_WARNING}The configuration path {'.'.join(keys[:-1])} is not present in original config. The variable {env_var_key} with value {env_var_value} can therefore not be set!{COLOR_END}", file=sys.stderr)

success = False
models = os.environ.get("CB_MODELS", "")
if models:
    models = models.split(",")

"""
user_machine = os.environ["CB_MACHINES"]

machine_models = {user_machine: 0}
to_remove = list()
for model_name, benchmark_config in config["benchmarks"].items():
    if models and model_name not in models:
        to_remove.append(model_name)
        print(f"{COLOR_WARNING}Skipping config generation for model {model_name} as requested by the user!{COLOR_END}", file=sys.stderr)
        continue

    for machine in list(benchmark_config["machine_parameters"].keys()):
        if machine != user_machine:
            del benchmark_config["machine_parameters"][machine]
        else:
            machine_models[machine] += 1
for model in to_remove:
    del config["benchmarks"][model]
"""
#for machine, count in machine_models.items():
#    if count == 0:
#        machines.remove(machine)
# Return list of machines from script to use it in the shell
#print(",".join(machines))

# Validate model parameters
for model_name, benchmark_config in config["benchmarks"].items():
    for param in benchmark_config.get("model_parameters", []):
        if "type" in param and not ("mode" in param and param["mode"] == "python"):
            try:
                if param["type"] == "float":
                    for c in param["_"].split(","):
                        float(c)
                elif param["type"] == "int":
                    for c in param["_"].split(","):
                        int(c)
            except ValueError:
                raise ValueError(f"Parameter {param['name']} for model {model_name} with value {param['_']} of type {param['type']} contains non-numeric values (e.g. variables), but 'mode' is not set to 'python'.")

# Write patched config file
with open(config_out_path, 'w') as config_outfile:
    yaml.dump(config, config_outfile)
