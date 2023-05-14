# Construct a jube-yaml file for a benchmark based on the configuration.
import sysy
import os

import yaml


jube_template_path = sys.argv[1]  # path to template of jube benchmark file
benchmark_config_path = sys.argv[2]  # path to configuration file of specific benchmark
benchmark_model_files_dir = sys.argv[3]  # path to directory of model files
benchmark_script_path = sys.argv[4]  # relative to benchmark_model_files_dir
benchmark_log_path = sys.argv[5]  # path where benchmark logs should be placed relative to workpackage dir

with(jube_template_path, 'r') as template_file:
    jube_benchmark = yaml.safe_load(template_file)

# Model parameters
with(benchmark_config_path, 'r') as benchmark_config_file:
    model_parameters = yaml.safe_load(benchmark_config_file)["model_parameters"]

# Name and outpath
benchmark_name = next(param["_"] for param in model_parameters if param["name"] == "model_name")
jube_benchmark["name"] = benchmark_name
jube_benchmark["outpath"] += benchmark_name

# Parameterset
jube_benchmark["parameterset"]["parameter"] = [
    {"name": "run_file", "type": "string", "_": "${jube_wp_abspath}" + os.sep + benchmark_script_path}
    {"name": "log_path", "type": "string", "_": "${jube_wp_abspath}" + os.sep + benchmark_log_path}
]

# Fileset
jube_benchmark["fileset"]["copy"] = ["$model_path" + os.sep + f for f in os.listdir(benchmark_model_files_dir) if os.path.isfile(os.path.join(benchmark_model_files_dir, f))]

# Substituteset
jube_benchmark["substituteset"]["iofile"] = {"in": benchmark_script_path, "out": benchmark_script_path}
jube_benchmark["substituteset"]["sub"] = list()
for param in model_parameters:
    jube_benchmark["substituteset"]["sub"].append({"source": '"{' + param["name"] + '}"', "dest": "$" + param["name"]})

# Bench step
next(u for u in next(b for b in jube_benchmark["step"] if b["name"] == "bench")["use"] if u["from"] == "MODEL_CONFIG")["from"] = benchmark_config_path

# Analyser and result
jube_benchmark["result"]["table"]["column"] = []
