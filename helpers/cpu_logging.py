import os
import sys
import json

save_path = sys.argv[1]

cpu_info = [element.strip().replace(' ', '').replace('(', '').replace(')', '')
            for element in os.popen('lscpu').readlines()]

cpu_info_dict = {}
for element in cpu_info:
    key, value = element.split(':')
    cpu_info_dict[key] = value

with open(os.path.join(save_path, 'cpu.json'), 'w') as f:
    json.dump(cpu_info_dict, f)
