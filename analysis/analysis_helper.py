import os
import yaml

import pickle

def shell(command):
    return os.system(command)


def shell_return(command):
    return os.popen(command).read().strip()


def load(filepath):
    with open(filepath, 'rb') as f:
        data = pickle.load(f)
    return data
    


def update_catalogue(catalogue_path, uuidgen_hash, cpu_info, job_info):
    dict_ = {
        uuidgen_hash: {
            'machine': os.popen('echo $HOSTNAME').read().strip(),
        }
    }

    dict_[uuidgen_hash].update(job_info)
    dict_[uuidgen_hash].update(cpu_info)

    with open(catalogue_path, 'r') as c:
        catalogue = yaml.safe_load(c)
        catalogue.update(dict_)

    with open(catalogue_path, 'w') as c:
        yaml.dump(catalogue, c)
        