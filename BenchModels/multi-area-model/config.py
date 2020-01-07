import os

script_dir = os.path.dirname(__file__)
multi_name = ''
simulation_name = 'simulations'
base_file_path = os.path.join(script_dir, multi_name)

# Absolute path of repository
base_path = base_file_path
print(base_path)
# Place to store simulations
data_path = os.path.join(base_file_path, simulation_name)
# Template for job scripts
jobscript_template = ''''''

# Command to submit jobs on the local cluster
submit_cmd = None

