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
yaml.add_representer(str, lambda dumper, data: dumper.represent_scalar('tag:yaml.org,2002:str', data, style='"'))

COLOR_WARNING = '\033[93m'
COLOR_FAIL = '\033[91m'
COLOR_END = '\033[0m'


jube_template_path = sys.argv[1]     # path to template of jube benchmark file
config_path = sys.argv[2]            # path to configuration file
builder_modules_dir = sys.argv[3]    # directory where builder installs modules
simulator_branch = sys.argv[4]       # repository branch of the simulator to benchmark
benchmark_models_dir = sys.argv[5]   # path to directory of model files for all models
benchmark_output_dir = sys.argv[6]   # path to directory where all jube benchmark files are placed
config_output_dir = sys.argv[7]      # path to directory where all benchmark config files are placed
scripts_dir = sys.argv[8]            # path to directory where all benchmark helper scripts are located
jube_results_dir = sys.argv[9]       # path to directory where all benchmark results should be placed
prepare_file_path = sys.argv[10]     # path to the prepare.sh file which is loaded on the compute node
pipeline_file_path = sys.argv[11]    # path to current pipeline definition file
machine = sys.argv[12]               # name of the machine the benchmark is running on

config_dir = os.path.dirname(config_path)

success = False

models = os.environ.get("CB_MODELS", "")
if models:
    models = models.split(",")

# Prepare output directory
if not os.path.exists(config_output_dir):
    os.makedirs(config_output_dir)

# Write prepare-file which is run right before the benchmark
with open(prepare_file_path, 'a') as prepare_file:
    with open(pipeline_file_path, 'r') as pipeline_file:
        pipeline = yaml.safe_load(pipeline_file)
    prepare_commands = pipeline.get(".environment-load", "")
    if prepare_commands:
        prepare_file.writelines([c.strip() + '\n' for c in prepare_commands])


with open(jube_template_path, 'r') as template_file:
    jube_benchmark_template = yaml.safe_load(template_file)
    jube_benchmark_template["outpath"] = jube_results_dir


with open(config_path, 'r') as config_file:
    config = yaml.safe_load(config_file)

for model_name, benchmark_config in config["benchmarks"].items():
    # Skip benchmarks defined in config which were not explicitly requested to be benchmarked
    if models and model_name not in models:
        print(f"{COLOR_WARNING}Skipping config file generation for model {model_name} as requested by the user!{COLOR_END}")
        continue

    if machine not in benchmark_config["machine_parameters"]:
        print(f"{COLOR_WARNING}Skipping config file generation for model {model_name} because no configuration was provided for machine {machine}!{COLOR_END}")
        continue
    
    jube_benchmark = copy.deepcopy(jube_benchmark_template)

    run_parameters = benchmark_config.get("run_parameters", dict())
    model_script_path = run_parameters["script_path"]
    metrics = benchmark_config.get("metrics", [])
    metrics_sum = benchmark_config.get("metrics_sum", [])
    all_metrics = metrics + metrics_sum

    # Write all metrics to a file to later extract timer data depending on the desired metrics
    metrics_out = {"metrics": metrics, "metrics_sum": metrics_sum}
    metrics_out_path = os.path.join(config_output_dir, f"cb_{model_name}_{machine}_metrics.yml")
    with open(metrics_out_path, 'w') as metrics_outfile:
        yaml.dump(metrics_out, metrics_outfile, default_flow_style=False)

    # Name
    jube_benchmark["name"] = model_name

    # Machine parameters
    machine_params = benchmark_config["machine_parameters"][machine]

    machine_param_list = list()
    for k in ["num_nodes", "tasks_per_node", "threads_per_task"]:
        for c in machine_params[k]:
            if not c.isdigit() or c != ",":
                machine_param_list.append({"name": k, "mode": "python", "type": "int", "_": machine_params[k]})
                break
        else:
            machine_param_list.append({"name": c, "type": "int", "_": machine_params[k]})
        w = machine_params["walltime"]
        if w[2] != ":" or w[5] != ":" or len(w) != 8 or not all(w[p].isdigit() for p in [0,1,3,4,6,7]):
            machine_param_list.append({"name": "walltime", "mode": "python", "type": "string", "_": w})
        else:
            machine_param_list.append({"name": "walltime", "type": "string", "_": w})
            
    jube_benchmark["parameterset"].append(
        {
            "name": "machine_parameters",
            "parameter": machine_param_list
        })
                

    jube_benchmark["parameterset"].append(
        {
            "name": "user_config",
            "parameter": [
                {"name": "account", "type": "string", "_": config["user_config"][machine]["account"]},
                {"name": "partition", "type": "string", "_": config["user_config"][machine]["partition"]}
            ]
        })

    # Software parameters
    jube_benchmark["parameterset"].append(
        {
            "name": "software_parameters",
            "parameter": [
                {"name": "simulator_branch", "_": simulator_branch},
                {"name": "buildermod", "_": builder_modules_dir}
            ]
        })

    # Model parameters
    jube_benchmark["parameterset"].append(
    {
        "name": "model_parameters",
        "parameter": benchmark_config.get("model_parameters", [])
    })

    # TODO: Logpath is currently hardcoded in models
    # Parameterset
    benchmark_model_files_dir = os.path.join(benchmark_models_dir, model_name)
    next(p for p in jube_benchmark["parameterset"] if p["name"] == "bench_parameters")["parameter"].extend([
        {"name": "run_file", "type": "string", "_": "${jube_wp_abspath}" + os.sep + model_script_path},
        {"name": "log_path", "type": "string", "_": "${jube_wp_abspath}" + os.sep + run_parameters.get("log_path", "")},
        {"name": "scripts_dir", "type": "string", "_": scripts_dir},
        {"name": "metrics_path", "type": "string", "_": metrics_out_path},
        {"name": "optional_run_command", "type": "string", "_": run_parameters.get("optional_run_command", "")},
        {"name": "run_args", "type": "string", "_": run_parameters.get("run_args", "")},
        {"name": "model_path", "type": "string", "_": os.path.abspath(benchmark_model_files_dir)},
        {"name": "prepare_file_path", "type": "string", "_": prepare_file_path}
    ])

    # Fileset
    jube_benchmark["fileset"]["copy"].extend([os.path.join(benchmark_model_files_dir, f) for f in os.listdir(benchmark_model_files_dir)])

    # Substituteset
    substitutions = next(s for s in jube_benchmark["substituteset"] if s["name"] == "bench_substitutions")
    substitutions["iofile"] = {"in": model_script_path, "out": model_script_path}
    for param in benchmark_config.get("model_parameters", []) + machine_params:
        substitutions["sub"].append({"source": "{" + param["name"] + "}", "dest": "$" + param["name"]})

    # Analyser
    jube_benchmark["patternset"] = {"name": "bench_patternset", "pattern": [
        {"name": key, "mode": "pattern", "type": "float", "_": key + " $jube_pat_fp"} for key in all_metrics
    ]}

    # Result
    jube_benchmark["result"]["table"]["column"] = benchmark_config.get("results", []) + all_metrics

    # Write jube benchmark file
    with open(os.path.join(benchmark_output_dir, f"cb_{model_name}_{machine}.yml"), 'w') as jube_outfile:
        yaml.dump(jube_benchmark, jube_outfile, default_flow_style=False)
    
    success = True

if not success:
    print(f"{COLOR_FAIL}No config for any model provided for this machine, aborting...{COLOR_END}", file=sys.stderr)
    exit(1)
