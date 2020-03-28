import os, re
from fnmatch import filter as fnfilter

# hook into snakemake logger
# be aware that this is not a default logging.Logger instance and calls behave
# a bit different.
from .logging import logger as log

# local workflow configuration
# edit the configfile to adjust the setup to your needs
configfile: "workflow.yaml"


rule install_jube:
    '''
    Install jube into the launch environment.

    This rule also creates a symlink to the launch environment, such that you
    can do
       $ conda activate environment/launch
    to enter the launch environment with your shell to check the correct setup.
    '''
    input:
        "JUBE-%s.tar.gz" % config['jubeversion'],
    output:
        info = ".jube-version-installed",
        env = directory('environment/launch'),
    conda: 'environment/launch.yaml'
    shell:
        '''
        set -x
        echo "CONDA_DEFAULT_ENV: $CONDA_DEFAULT_ENV"
        echo "CONDA_PREFIX: $CONDA_PREFIX"
        echo "installing to CONDA_PREFIX..."
        pip install --prefix $CONDA_PREFIX {input}   # this uses the conda pip to install into the correct environment
        jube --version >{output} 2>&1
        echo "installed in $CONDA_PREFIX" >>{output.info}
        ln -vs $CONDA_PREFIX {output.env}
        '''

rule download_jube:
    '''
    (called from 'install_jube')
    Fetch a given JUBE version from the source website.
    '''
    output:
        temporary("JUBE-{version}.tar.gz")
    shell:
        '''
        wget http://apps.fz-juelich.de/jsc/jube/jube2/download.php?version={wildcards.version} -O JUBE-{wildcards.version}.tar.gz
        '''

def bench_ids():
    # construct the list of available bench_runs
    # (this is required e.g. for rule show_runs)
    log.info("searching for available runs in %s" % config['outpath'])
    UUID_re = r'[0-9a-f]{8}-([0-9a-f]{4}-){3}[0-9a-f]{12}'
    ids = [name for name in os.listdir(config['outpath']) if re.match(UUID_re, name) and os.path.isdir(os.path.join(config['outpath'], name))]
    log.info("%d bench_runs found." % len(ids))
    return ids

rule show_help:
    '''
    (default action)
    Print some usage information and give pointers to where to find more
    documentation.
    '''
    shell:
        '''read documentation << :


Benchmarking Workflow
=====================

The overall workflow is configured in the config file used for snakemake. This
is usually the 'workflow.yaml' or a similar file given with --config option to
Snakemake. Look into that file to see a description of the available config
options.

Once this is configured to your needs, several actions can be performed in the
workflow. Run each action simply with 'snakemake <action>', e.g. 'snakemake run'.
To see a full list of available actions and a short description run:

    snakemake --list

To look at a graph of triggered actions, e.g. for the 'run' command:

    snakemake --dag run | dot -Tx11

Or in case of many nodes being involved 

    snakemake --rulegraph show_runs | dot -Tx11


For more in-depth documentation and how to add new models / machines /
run-types etc. read the corresponding files in doc/.

:
false # intentional stop
        '''

rule run:
    '''
    create a new bench_run according to the default setting in workflow.yaml
    '''
    message:
        '''

################################################################################
                           R U N   C O M P L E T E
################################################################################

        '''
    input:
        '{model}__{benchmark}'.format(**config['default'])
    shell:
        '''
        echo "MOVING LOGFILES INTO THE RUN DIRECTORY"
        mkdir -v $(cut -f1 -d' ' {input}.submit.id)/logs
        mv -v *.log $(cut -f1 -d' ' {input}.submit.id)/logs
        echo "REMOVING ID FILE OF COMPLETED RUN:"
        cat {input}.submit.id
        rm -f {input}.submit.id
        '''


rule cancel:
    '''
    remove all local files from a broken run (does not stop queued jobs)
    '''
    shell:
        '''
        for runfile in $(ls -1 *.submit.id 2>/dev/null || echo "No run found to cancel!" >&2); do
            echo "WILL REMOVE $runfile"
            echo "WILL REMOVE"
            du -hcs "$(cat $runfile | cut -f1 -d' ')"/* || true;
            echo "Sure? (type 'yes')?"
            read sure
            test "$sure" == "yes"
            rm -rf "$(cat $runfile | cut -f1 -d' ')"
            rm -f "$runfile"
        done
        '''

rule info:
    '''
    show jube info of current run
    '''
    input:
        runs = lambda wc: fnfilter(os.listdir(), "*.submit.id"),
        jube = ".jube-version-installed",
    conda:
        'environment/launch.yaml'
    shell:
        '''
        if [ -z "{input.runs}" ]; then
            echo "No current runs found."
            exit 0;
        fi
        for ID in {input.runs}; do
            echo ""
            echo ""
            echo "$ID"
            jube info $(cut -f1 -d' ' "$ID")
            jube info $(cat "$ID")
        done
        '''
rule submit_run:
    '''
    (called from 'run')
    Start a new benchmark run and hand control over to JUBE.

    In case something goes wrong here, you can clean-up the broken run by

      snakemake cancel

    '''
    input:
        jubefile = 'test_bench.xml',
        jube = ".jube-version-installed",
    output:
        protected('{model}__{benchmark,[^.]+}.submit.id'),
    log:
        '{model}__{benchmark,[^.]+}.submit.log',
    conda:
        'environment/launch.yaml'
    message:
        '''

################################################################################
                S T A R T I N G   N E W   J U B E   R U N
################################################################################

        '''
    shell:
        '''
        jube -v run {input.jubefile} --include-path machine/{config[default][machine]} benchmark/{config[default][benchmark]} model/{config[default][model]} --outpath {config[outpath]}/$(uuidgen) |& tee {log}
        common/submit-info.py {log} >{output}
        '''


rule continue_run:
    '''
    (called from 'run')
    Tell JUBE to continue the last started benchmark run.
    '''
    input:
        id = '{model}__{benchmark}.submit.id',
    output:
        temporary('{model}__{benchmark}.complete'),
    log:
        run = '{model}__{benchmark,[^.]+}.continue.log',
        id = '{model}__{benchmark,[^.]+}.info.log',
    conda:
        'environment/launch.yaml'
    shell:
        '''
        set -x
        jube info $(cat {input.id}) >{log.id} 2>&1
        date >>{log.run}
        jube continue $(cat {input.id}) |& tee -a {log.run}
        if grep "Benchmark finished" {log.id}; then
            cp {input.id} {output}
        else
            echo -e "\n\n BENCHMARK NOT FINISHED! (rerun "smake run" to continue)\n\n"
        fi
        '''

rule result_aggregation:
    '''
    (called from 'run')

    Tell JUBE to fetch 'results' and run 'analyse'.
    '''
    input:
        id = '{model}__{benchmark}.complete',
    output:
        temporary('{model}__{benchmark,[^.]+}.analysed'),
    log:
        result = '{model}__{benchmark,[^.]+}.result.log',
    message:
        '''
        Analysing the runs and creating result files...

        This may take a while.

        Check analyse.log in output directory.
        '''
    conda:
        'environment/launch.yaml'
    shell:
        '''
        jube analyse $(cat {input.id})
        ID=$(cut -f1 -d' ' {input.id})
        mkdir -pv $ID/results/txt $ID/results/csv
        egrep -o '<table.*' test_bench.xml | egrep -o 'name="[^"]+"' | egrep -o '[^"]+' | grep -v name |\
        while read tablename; do
            echo "Result table $tablename..."
            jube result $ID -o $tablename --style pretty >$ID/results/txt/$tablename
            jube result $ID -o $tablename --style csv    >$ID/results/csv/$tablename
        done
        echo "Results stored in $ID/results/"
        cp {input.id} {output}
        '''


rule result_cleanup:
    '''
    (called from 'run')
    Run `JUBE analyse` and `JUBE result` to compile result dataset. Then tidy up.
    '''
    input:
        '{model}__{benchmark}.analysed',
    output:
        temporary('{model}__{benchmark,[^.]+}.tidy'),
    shell:
        '''
        UUID=$(cat {input})
        echo "Annotating $UUID"
        echo git annex metadata model={wildcards.model} $UUID
        echo git annex metadata benchmark={wildcards.benchmark} $UUID
        echo git annex metadata model={wildcards.model} $UUID
        echo $UUID >{output}
        '''

rule upload_results:
    '''
    (called from 'run')
    Add metadata to tidy results, add them to git annex and sync out to the
    central repository.
    '''
    input:
        '{model}__{benchmark}.tidy',
    output:
        temporary('{model}__{benchmark,[^.]+}'),
    message:
        '''
        ALL RESULTS FOR {output} ARE PREPARED AND READY TO BE UPLOADED.
        '''
    shell:
        '''
        UUID=$(cat {input})
        cat {input}
        echo git annex add $UUID
        echo git ci -m "added {wildcards.benchmark} of {wildcards.model} on {config[default][machine]}"
        touch {output}
        '''

rule show_runs:
    '''
    Create overviews of all available outputs and give very brief status of
    each.
    '''
    input:
        os.path.join(config['outpath'], "_OVERVIEW_.txt")
    shell:
        '''
        cat {input}
        '''

rule run_overview:
    '''
    (called from 'run_overview')
    Creates a summary by extracting lines from extract_run_info results.
    '''
    input:
        expand("{path}/{uuid}.info.txt", path=config['outpath'], uuid=bench_ids()),
    output:
        os.path.join(config['outpath'], "_OVERVIEW_.txt")
    shell:
        '''
        #grep -H "Benchmark.*finished" {input}
        egrep -H "^ +[0-9]+|Benchmark.*finished" output/*.info.txt | sed -e "s/^.*--- Benchmark/|/g" -e "s/ ---$/%/" | tr "\n" " " | tr "%" "\n"  | sed -e "s/\.info\.txt:/ |/" -e "s%^.*/%%;" >{output}
        '''

rule extract_run_info:
    '''
    (called from 'run_overview')
    Create an info file for a specific bench_run (UUID)
    '''
    output:
        '{uuid}.info.txt',
    conda:
        'environment/launch.yaml'
    shell:
        '''
        jube info {wildcards.uuid} >{output}
        jube info {wildcards.uuid} --id 0 >>{output}
        '''
    
