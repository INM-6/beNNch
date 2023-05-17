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
Create a single configuration file for a specific benchmark by collecting all required information from various places.
"""

import sys
import os
import copy
import subprocess
import shutil

import yaml


config_file_path = sys.argv[1]
config_output_dir = sys.argv[2]
builder_modules_dir = sys.argv[3]
simulator_branch = sys.argv[4]
machine = sys.argv[5]

config_dir = os.path.dirname(config_file_path)


# Get all configurations which should be overwritten from environment variables
config_overwrites = list()
for env_var_key, env_var_value in os.environ.items():
    if env_var_key.startswith("CB_"):
        # TODO: Think about a good format for overwrites
        config_overwrites.append((env_var_key[3:].split("__"), env_var_value))


with open(config_file_path) as config_file:
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
        print("Skipping config file generation for model " + model_name)
        continue

    combined_configs[model_name] = dict()
    combined_config = {"parameterset": list()}

    # Software parameters
    combined_config["parameterset"].append(
        {
            "name": "software_parameters",
            "parameter": [
                {"name": "simulator", "_": "nest-simulator"},
                {"name": "version", "_": simulator_branch},
                {"name": "variant", "_": "default"},
                {"name": "suffix", "_": ""},
                {"name": "buildermod", "_": builder_modules_dir}
            ]
        })

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
            print("Skipping config file generation for machine " + machine_name)
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

        # User config
        combined_config["parameterset"].append(
            {
                "name": "user_config",
                "parameter": [
                    {"name": "account", "type": "string", "_": config["user_config"][machine]["account"]},
                    {"name": "email_address", "type": "string", "_": config["user_config"][machine]["email_address"]},
                    {"name": "partition", "type": "string", "_": config["user_config"][machine]["partition"]}
                ]
            })
        
        combined_config["metrics"] = model_config["metrics"]
        combined_config["metrics_sum"] = model_config["metrics_sum"]

        # TODO: Patch configuration using environment variables
        # to_remove = list()
        # for i, overwrite in enumerate(config_overwrites):
        #     try:
        #         to_overwrite = combined_config["parameterset"]
        #         for key_part in overwrite[0][:-1]:
        #             to_overwrite = to_overwrite[key_part]
        #         to_overwrite[key_part[-1]] = overwrite[1]
        #         to_remove.append(i)
        #     except:
        #         pass
        # # Remove all overwrites which have been successfully applied
        # for i in to_remove:
        #     del config_overwrites[i]

        # Dump config file
        with open(os.path.join(config_output_dir, os.path.basename(model_config_path)), 'w') as outfile:
            yaml.dump(combined_config, outfile, default_flow_style=False)

for conf in combined_configs.values():
    if len(conf) > 0:
        break
else:
    print("No config for any model provided for this machine, aborting...")
    exit(1)

# Check for any leftover (unapplied) overwrites
for overwrite_key, overwrite_value in config_overwrites:
    print(f"Benchmark variable {overwrite_key} with value {overwrite_value} unused and will be ignored!", file=sys.stderr)
