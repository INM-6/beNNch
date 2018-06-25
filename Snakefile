import os

configfile: "config.yaml"

rule all:
    input:
        '.JUBE.installed'

rule install_jube:
    input:
        "JUBE-%s.tar.gz" % config['jubeversion'],
    output:
        ".JUBE.installed"
    params:
        prefix = os.environ.get("CONDA_PREFIX", os.path.expanduser("~/.local"))
    shell:
        '''
        tar -xf {input}
        echo "installing to environment '$CONDA_DEFAULT_ENV'..."
        pip install ./JUBE-{config[jubeversion]}   # this uses the conda pip to install into the correct environment
        jube --version >{output}
        rm -rf JUBE-{config[jubeversion]}
        '''

rule download_jube:
    output:
        temporary("JUBE-{version}.tar.gz")
    shell:
        '''
        wget http://apps.fz-juelich.de/jsc/jube/jube2/download.php?version={wildcards.version} -O JUBE-{wildcards.version}.tar.gz
        '''
