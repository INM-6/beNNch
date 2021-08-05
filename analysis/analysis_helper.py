import os
import json
import yaml

import pickle


def shell(command):
    return os.system(command + '>/dev/null 2>&1')


def shell_return(command):
    return os.popen(command).read().strip()


def load(filepath):
    with open(filepath, 'rb') as f:
        data = pickle.load(f)
    return data

def git_annex(cpu_info, job_info, uuidgen_hash, base_path, result_path):

    tmp_result_file_path = os.path.join(base_path, uuidgen_hash + '.csv')
    result_file_path = os.path.join(result_path, uuidgen_hash + '.csv')

    machine = os.popen('echo $HOSTNAME').read().strip()

    shell(f'cp {tmp_result_file_path} {result_file_path}')
    shell(f'git annex add {result_file_path}')
    shell(f'git annex metadata {result_file_path} --set key={uuidgen_hash}')

    for info_dict in [job_info, cpu_info]:
        for key, value in info_dict.items():
            shell(f'git annex metadata {result_file_path} '
                  + f'--set {key}="{value}" --force')
    shell(f'git annex metadata {result_file_path} '
          + f'--set machine="{machine}" --force')


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
