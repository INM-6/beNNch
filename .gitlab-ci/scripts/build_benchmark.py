# Construct a jube-yaml file for a benchmark based on the configuration.
import sysy
import os

import yaml


jube_template_path = sys.argv[1]
benchmark_config_path = sys.argv[2]
benchmark_script_path = sys.argv[3]
benchmark_log_path = sys.argv[4]

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

# Substituteset

# Bench step

# Analyser and result
