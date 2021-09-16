import os
import sys
import pickle

save_path = sys.argv[1]
vps_per_node = sys.argv[2]
num_nodes = sys.argv[3]
tasks_per_node = sys.argv[4]
threads_per_task = sys.argv[5]
nest = sys.argv[6]
model_name = sys.argv[7]
network_state = sys.argv[8]
record_spikes = sys.argv[9]
affinity = ';'.join(sys.argv[10:])
job_info = {
    'vps_per_node': vps_per_node,
    'num_nodes': num_nodes,
    'tasks_per_node': tasks_per_node,
    'threads_per_task': threads_per_task,
    'affinity': affinity,
    'nest': nest,
    'model_name': model_name,
    'network_state': network_state,
    'record_spikes': record_spikes
}

with open(os.path.join(save_path, 'job.pkl'), 'wb') as f:
    pickle.dump(job_info, f)
