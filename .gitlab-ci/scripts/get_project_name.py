import sys
import os

import yaml


config_dir = sys.argv[1]
machine = sys.argv[2]


with open(os.path.join(config_dir, 'config.yml')) as config_file:
    config = yaml.safe_load(config_file)

project_name = config["user_config"][machine]["project"]
runner_account = config["user_config"][machine]["account"]

# Check if project name should be overwritten from environment variable
config_overwrites = list()
for env_var_key, env_var_value in os.environ.items():
    if env_var_key == f"CB_user_config__{machine}__project":
        project_name = env_var_value
    if env_var_key == f"CB_user_config__{machine}__account":
        runner_account = env_var_value

os.environ["PROJECT_NAME"] = project_name
os.environ["RUNNER_ACCOUNT"] = runner_account
