import os
import sys
import pickle

save_path = sys.argv[1]
base_vp = sys.argv[2]
num_vps = sys.argv[3]
num_nodes = sys.argv[4]
tasks_per_node = sys.argv[5]
threads_per_task = sys.argv[6]
nest = sys.argv[7]
network_state = sys.argv[8]

job_info = {
    'base_vp': base_vp,
    'num_vps': num_vps,
    'num_nodes': num_nodes,
    'tasks_per_node': tasks_per_node,
    'threads_per_task': threads_per_task,
    'nest': nest,
    'network_state': network_state,
}

with open(os.path.join(save_path, 'job.pkl'), 'wb') as f:
    pickle.dump(job_info, f)
