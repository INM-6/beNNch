import sys
import os

import yaml


requirement_type = sys.argv[1]
config_file_path = sys.argv[2]
machine = sys.argv[3]

with open(config_file_path) as config_file:
    config = yaml.safe_load(config_file)

requirements = config["requirements"][machine].get(requirement_type, [])

print(*requirements, sep=" ", end="")
