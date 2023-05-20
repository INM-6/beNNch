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
import copy

import yaml


jube_template_path = sys.argv[1]     # path to template of jube benchmark file
config_path = sys.argv[2]            # path to configuration file of specific benchmark
benchmark_models_dir = sys.argv[3]   # path to directory of model files for all models
benchmark_output_dir = sys.argv[4]   # path to directory where all jube benchmark files are placed
config_output_dir = sys.argv[5]      # path to directory where all benchmark config files are placed
scripts_dir = sys.argv[6]            # path to directory where all benchmark helper scripts are located
jube_results_dir = sys.argv[7]       # path to directory where all benchmark results should be placed
machine = sys.argv[8]                # name of the machine the benchmark is running on

config_dir = os.path.dirname(config_path)

with open(jube_template_path, 'r') as template_file:
    jube_benchmark_template = yaml.safe_load(template_file)
    jube_benchmark_template["outpath"] = jube_results_dir

models = os.environ.get("CB_MODELS", [])
if models:
    models = models.split(",")

with open(config_path, 'r') as config_file:
    config = yaml.safe_load(config_file)["benchmarks"]

for model_name, benchmark_configs in config.items():
    if models and model_name not in models:
        print("Skipping benchmark file generation for model " + model_name)
        continue

    model_config_path = os.path.join(config_dir, benchmark_configs["config_path"])

    # Model parameters
    with open(model_config_path, 'r') as benchmark_config_file:
        benchmark_config = yaml.safe_load(benchmark_config_file)
        model_parameters = benchmark_config["model_parameters"]
        run_parameters = benchmark_config["run_parameters"]
        metrics = benchmark_config["metrics"]
        metrics_sum = benchmark_config["metrics_sum"]

    model_script_path = run_parameters["script_path"]

    jube_benchmark = copy.deepcopy(jube_benchmark_template)

    # Name and outpath
    jube_benchmark["name"] = model_name

    # TODO: Logpath is currently hardcoded in models
    # Parameterset
    benchmark_model_files_dir = os.path.join(benchmark_models_dir, model_name)
    next(p for p in jube_benchmark["parameterset"] if p["name"] == "bench_parameters")["parameter"].extend([
        {"name": "run_file", "type": "string", "_": "${jube_wp_abspath}" + os.sep + model_script_path},
        {"name": "log_path", "type": "string", "_": "${jube_wp_abspath}" + os.sep + run_parameters.get("log_path", "")},
        {"name": "scripts_dir", "type": "string", "_": scripts_dir},
        {"name": "optional_run_command", "type": "string", "_": run_parameters.get("optional_run_command", "")},
        {"name": "run_args", "type": "string", "_": run_parameters.get("run_args", "")},
        {"name": "model_path", "type": "string", "_": os.path.abspath(benchmark_model_files_dir)},
        {"name": "model_config_path", "type": "string", "_": os.path.join(config_output_dir, os.path.basename(model_config_path))}
    ])

    # Fileset
    jube_benchmark["fileset"]["copy"].extend([os.path.join(benchmark_model_files_dir, f) for f in os.listdir(benchmark_model_files_dir) if os.path.isfile(os.path.join(benchmark_model_files_dir, f))])

    # Substituteset
    substitutions = {"name": "bench_substitutions", "iofile": {"in": model_script_path, "out": model_script_path}, "sub": list()}
    for param in model_parameters:
        substitutions["sub"].append({"source": "{" + param["name"] + "}", "dest": "$" + param["name"]})
    # TODO: How to deal with scale properly?
    substitutions["sub"].append({"source": "{scale}", "dest": "$scale"})
    substitutions["sub"].append({"source": "{num_vps}", "dest": "$num_vps"})
    jube_benchmark["substituteset"].append(substitutions)

    # Bench step
    step = next(b for b in jube_benchmark["step"] if b["name"] == "bench")
    next(u for u in step["use"] if type(u) == dict and u["from"] == "MODEL_CONFIG")["from"] = os.path.join(config_output_dir, os.path.basename(model_config_path))

    # Analyser and result
    jube_benchmark["result"]["table"]["column"] = metrics + metrics_sum

    with open(os.path.join(benchmark_output_dir, f"cb_{model_name}_{machine}.yml"), 'w') as outfile:
        yaml.dump(jube_benchmark, outfile, default_flow_style=False) 
