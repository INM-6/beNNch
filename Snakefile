import os

configfile: "config.yaml"

rule all:
    input:
        '.jube-version-installed'

rule install_jube:
    input:
        "JUBE-%s.tar.gz" % config['jubeversion'],
    output:
        ".jube-version-installed"
    conda: 'environment.yml'
    shell:
        '''
        echo "CONDA_DEFAULT_ENV: $CONDA_DEFAULT_ENV"
        echo "CONDA_PREFIX: $CONDA_PREFIX"
        echo "installing to CONDA_PREFIX..."
        pip install --prefix $CONDA_PREFIX {input}   # this uses the conda pip to install into the correct environment
        jube --version >{output} 2>&1
        echo "installed in $CONDA_PREFIX" >>{output}
        '''

rule download_jube:
    output:
        temporary("JUBE-{version}.tar.gz")
    shell:
        '''
        wget http://apps.fz-juelich.de/jsc/jube/jube2/download.php?version={wildcards.version} -O JUBE-{wildcards.version}.tar.gz
        '''
