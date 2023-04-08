import sys
import os
import copy
import subprocess
import shutil

import yaml


config_dir = sys.argv[1]
jube_output_dir = sys.argv[2]
config_output_dir = sys.argv[3]
builder_modules_dir = sys.argv[4]
bennch_dir = sys.argv[5]
model_repos_dir = sys.argv[6]
machine = sys.argv[7]


# Get all configurations which should be overwritten from environment variables
config_overwrites = list()
for env_var_key, env_var_value in os.environ.items():
    if env_var_key.startswith("CB_"):
        # TODO: Think about a good format for overwrites
        config_overwrites.append((env_var_key[4:].split("__"), env_var_value))


with open(os.path.join(config_dir, 'config.yml')) as config_file:
    config = yaml.safe_load(config_file)

combined_configs = dict()

models = os.environ.get("CB_MODELS", [])
if models:
    models = models.split(",")

# Prepare output directory
if not os.path.exists(config_output_dir):
    os.makedirs(config_output_dir)

for model_name, benchmark_configs in config["benchmarks"].items():
    # Skip benchmarks defined in config which were not explicitly requested to be benchmarked
    if models and model_name not in models:
        continue

    combined_configs[model_name] = dict()

    # Clone model repository into temporary location
    model_temp_dir = os.path.join(model_repos_dir, "tmp", model_name)
    subprocess.run(["git", "clone", "--depth=1", "--single-branch", "--branch=" + benchmark_configs["branch"], benchmark_configs["repository"], model_temp_dir])
    # Copy relevant part of the repository to final destination
    model_repo_dir = os.path.join(model_repos_dir, model_name)
    shutil.copytree(os.path.join(model_temp_dir, benchmark_configs["model_path"]), model_repo_dir)
    # Remove temporary repository again
    shutil.rmtree(model_temp_dir)

    # File paths
    combined_config = {"parameterset": [
    {
        "name": "file_paths",
        "parameter": [
            {"name": "base_path", "type": "string", "_": os.path.abspath(bennch_dir)},
            {"name": "model_path", "type": "string", "_": os.path.abspath(model_repo_dir)}
        ]
    }]}

    # Software parameters
    combined_config["parameterset"].append(
        {
            "name": "software_parameters",
            "parameter": [
                {"name": "simulator", "_": "nest-simulator"},
                {"name": "version", "_": "benchmarking"},
                {"name": "variant", "_": "default"},
                {"name": "suffix", "_": ""},
                {"name": "buildermod", "_": builder_modules_dir}
            ]
        })
    
    model_jube_path = os.path.join(config_dir, benchmark_configs["jube_path"])
    model_config_path = os.path.join(config_dir, benchmark_configs["config_path"])
    with open(model_config_path) as model_config_file:
        model_config = yaml.safe_load(model_config_file)

    # Model parameters
    combined_config["parameterset"].append(
    {
        "name": "model_parameters",
        "parameter": model_config["model_parameters"]
    })

    for machine_name, machine_config_path in model_config["machines"].items():
        # Only generate config for the machine this script runs on
        if machine_name != machine:
            continue

        combined_configs[model_name][machine_name] = copy.deepcopy(combined_config)

        with open(os.path.join(os.path.dirname(model_config_path), machine_config_path)) as machine_config_file:
            machine_config = yaml.safe_load(machine_config_file)

        # Machine parameters
        combined_config["parameterset"].append(
            {
                "name": "machine_parameters",
                "parameter": machine_config["machine_parameters"]
            })

        # Patch configuration using environment variables
        to_remove = list()
        for i, overwrite in enumerate(config_overwrites):
            try:
                to_overwrite = combined_config["parameterset"]
                for key_part in overwrite[0][:-1]:
                    to_overwrite = to_overwrite[key_part]
                to_overwrite[key_part[-1]] = overwrite[1]
                to_remove.append(i)
            except:
                pass
        # Remove all overwrites which have been successfully applied
        for i in to_remove:
            del config_overwrites[i]

        # Write out configuration and jube files
        # Dump Jube file
        if not os.path.exists(jube_output_dir):
            os.makedirs(jube_output_dir)
        shutil.copyfile(model_jube_path, os.path.join(jube_output_dir, f"cb_{model_name}_{machine_name}.yml"))
        # Dump config file
        with open(os.path.join(config_output_dir, os.path.basename(model_config_path)), 'w') as outfile:
            yaml.dump(combined_config, outfile, default_flow_style=False)

# Generate, patch, and dump user config
user_config = {"parameterset": {"name": "user_config", "parameter": 
    [{"name": "account", "type": "string", "_": config["user_config"][machine]["account"]},
    {"name": "email_address", "type": "string", "_": config["user_config"][machine]["email_address"]},
    {"name": "partition", "type": "string", "_": config["user_config"][machine]["partition"]}]
}}
    

# Patch configuration using environment variables
to_remove = list()
for i, overwrite in enumerate(config_overwrites):
    try:
        to_overwrite = user_config
        for key_part in overwrite[0][:-1]:
            to_overwrite = to_overwrite[key_part]
        to_overwrite[key_part[-1]] = overwrite[1]
        to_remove.append(i)
    except:
        pass
for i in to_remove:
    del config_overwrites[i]

with open(os.path.join(config_output_dir, "user_config.yaml"), 'w') as outfile:
    yaml.dump(user_config, outfile, default_flow_style=False) 

# Check for any leftover (unapplied) overwrites
for overwrite_key, overwrite_value in config_overwrites:
    print(f"Benchmark variable {overwrite_key} with value {overwrite_value} unused and will be ignored!", file=sys.stderr)
