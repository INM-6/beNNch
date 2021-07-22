import os
import yaml

def shell(command):
    return os.system(command)


def shell_return(command):
    return os.popen(command).read().strip()


def update_catalogue(catalogue_path, uuidgen_hash):
    dict_ = {
        uuidgen_hash: {
            'machine': os.popen('echo $HOSTNAME').read().strip(),
            'cpu model name': shell_return('grep -m 1 "model name" /proc/cpuinfo'),
            'notes': [
                {'num vps per node': 128},
                {'MPI process per node': 4},
                {'threads per MPI proc': 32},
                {'nest': '2.14.1 with timers'},
                {'mam_state': 'Fig5'}
            ],
            'plot_name': 'scaling_2_14_1_Fig5',
            'reason': 'Benchmark comparison with NEST 3.',
            'where': [
                'jureca.fz-juelich.de',
                '/p/project/cjinb33/jinb3330/gitordner/BenchWork/jube_MAM/nest_2/000003'
            ]
        }
    }

    with open(catalogue_path, 'r') as c:
        catalogue = yaml.safe_load(c)
        catalogue.update(dict_)

    with open(catalogue_path, 'w') as c:
        yaml.dump(catalogue, c)

    # if overwrite:
    #     catalogue_fn = 'catalogue.yaml'
    #     with open('catalogue.yaml', 'w') as yamlfile:
    #         yaml.safe_dump(catalogue, yamlfile)
    # else:
    #     catalogue_fn = 'catalogue_new.yaml'
    #     with open('catalogue_new.yaml', 'w') as yamlfile:
    #         yaml.safe_dump(catalogue, yamlfile)