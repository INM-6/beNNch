import os
import json

import numpy as np
import pandas as pd


def shell(command):
    return os.system(command)


def shell_without_print(command):
    return os.system(command + '>/dev/null 2>&1')


def shell_return(command):
    return os.popen(command).read().strip()


def load(filepath):
    with open(filepath, 'r') as f:
        data = json.load(f)
    return data


def git_annex(cpu_info, job_info, uuidgen_hash, base_path):

    tmp_result_file_path = os.path.join(base_path, uuidgen_hash + '.csv')
    result_file_path = os.path.join('./', uuidgen_hash + '.csv')
    tmp_archive_path = os.path.join(base_path, 'metadata_archive')
    archive_path = os.path.join('./', f'{uuidgen_hash}.tar.gz')

    shell(f'tar -czf {archive_path} -C {tmp_archive_path} .')

    # works for machines with the naming scheme XXX.name (used for JSC
    # clusters, might need adjustment for other machines)
    machine = os.popen('echo $HOSTNAME').read().strip().split('.')[-1]
    user = os.popen('echo $USER').read().strip()

    shell(f'cp {tmp_result_file_path} {result_file_path}')
    shell(f'git annex add {result_file_path}')
    shell(f'git annex add {archive_path}')
    shell_without_print(
        f'git annex metadata {result_file_path} --set key={uuidgen_hash}')

    for info_dict in [job_info, cpu_info]:
        for key, value in info_dict.items():
            value_without_spaces = value.replace(' ', ';')
            shell_without_print(f'git annex metadata {result_file_path} '
                                + f'--set {key}="{value_without_spaces}" --force')
    # additionally add machine and user name
    shell_without_print(f'git annex metadata {result_file_path} '
                        + f'--set machine="{machine}" --force')
    shell_without_print(f'git annex metadata {result_file_path} '
                        + f'--set user="{user}" --force')

    averaged_over = len(
        np.unique(pd.read_csv(result_file_path)['rng_seed'].values))
    shell_without_print(f'git annex metadata {result_file_path} '
                        + f'--set averaged_over="{averaged_over}" --force')
