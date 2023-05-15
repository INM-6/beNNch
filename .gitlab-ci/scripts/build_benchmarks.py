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
Construct a jube-yaml file for a benchmark based on the configuration.
"""

import sys
import os

import yaml


jube_template_path = sys.argv[1]     # path to template of jube benchmark file
config_path = sys.argv[2]            # path to configuration file of specific benchmark
benchmark_models_dir = sys.argv[3]   # path to directory of model files for all models
benchmark_output_dir = sys.argv[4]   # path to directory where all jube benchmark files are placed
scripts_dir = sys.argv[5]            # path to directory where all benchmark helper scripts are located
machine = sys.argv[6]                # name of the machine the benchmark is running on

config_dir = os.path.dirname(config_file_path)

with(jube_template_path, 'r') as template_file:
    jube_benchmark_template = yaml.safe_load(template_file)

models = os.environ.get("CB_MODELS", [])
if models:
    models = models.split(",")

# Prepare output directory
if not os.path.exists(config_output_dir):
    os.makedirs(config_output_dir)

with(config_path, 'r') as config_file:
    config = yaml.safe_load(config_file)["benchmarks"]

for model_name, benchmark_configs in config.items():
    if models and model_name not in models:
        print("Skipping benchmark file generation for model " + model_name)
        continue

    model_config_path = os.path.join(config_dir, benchmark_configs["config_path"])

    # Model parameters
    with(model_config_path, 'r') as benchmark_config_file:
        benchmark_config = yaml.safe_load(benchmark_config_file)
        model_parameters = benchmark_config["model_parameters"]
        run_parameters = benchmark_config["run_parameters"]
        metrics = benchmark_config["metrics"]
        metrics_sum = benchmark_config["metrics_sum"]

    model_script_path = run_parameters["script_path"]

    # Name and outpath
    jube_benchmark["name"] = next(param["_"] for param in model_parameters if param["name"] == "model_name")

    # TODO: Logpath is currently hardcoded in models
    # Parameterset
    next(p for p in jube_benchmark["parameterset"] if p["name"] == "bench_parameters")["parameter"].extend([
        {"name": "run_file", "type": "string", "_": "${jube_wp_abspath}" + os.sep + model_script_path},
        {"name": "log_path", "type": "string", "_": "${jube_wp_abspath}" + os.sep + run_parameters.get("log_path", "")},
        {"name": "scripts_dir", "type": "string", "_": scripts_dir},
        {"name": "optional_run_command", "type": "string", "_": run_parameters.get("optional_run_command", "")},
        {"name": "run_args", "type": "string", "_": run_parameters.get("run_args", "")},
        {"name": "model_path", "type": "string", "_": os.path.abspath(benchmark_model_files_dir)}
    ])

    # Fileset
    jube_benchmark["fileset"]["copy"].extend(["$model_path" + os.sep + f for f in os.listdir(benchmark_model_files_dir) if os.path.isfile(os.path.join(benchmark_model_files_dir, f))])

    # Substituteset
    substitutions = {"name": "bench_substitutions", "iofile": {"in": model_script_path, "out": model_script_path}, "sub": list()}
    for param in model_parameters:
        substitutions["sub"].append({"source": '"{' + param["name"] + '}"', "dest": "$" + param["name"]})
    jube_benchmark["substituteset"].append(substitutions)

    # Bench step
    next(u for u in next(b for b in jube_benchmark["step"] if b["name"] == "bench")["use"] if u["from"] == "MODEL_CONFIG")["from"] = benchmark_config_path

    # Analyser and result
    jube_benchmark["result"]["table"]["column"] = metrics + metrics_sum

    with open(os.path.join(benchmark_output_dir, f"cb_{model_name}_{machine_name}.yml"), 'w') as outfile:
        yaml.dump(user_config, outfile, default_flow_style=False) 
