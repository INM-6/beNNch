

# conda-environment-stacks

    conda env create -f env1.yaml -p ./env1
    conda env create -f env2.yaml -p ./env2

   When entering one of the environments, you have the corresponding tools available:

    (base) $ conda activate ./env1

    (env1) $ snakemake
    5.14.0

    (env1) $ dot -V
    -bash: dot: command not found

   In the other environment vice versa

    (base) $ conda activate ./env2

    (env2) $ snakemake
    -bash: snakemake: command not found

    (env2) $ dot -V
    dot - graphviz version 2.40.1 (20161225.0304)

   The interesting part is loading one env after another:

    (base) $ conda activate ./env1
    (env1) $ conda activate ./env2

    (env2) $ snakemake
    -bash: snakemake: command not found

    (env2) $ dot -V
    dot - graphviz version 2.40.1 (20161225.0304)

   while

    (base) $ conda activate ./env1
    (env1) $ conda activate --stack ./env2

    (env2) $ snakemake
    5.14.0

    (env2) $ dot -V
    dot - graphviz version 2.40.1 (20161225.0304)

   Interesting things to look at are:

   * `env | grep "CONDA_"`
   * `which snakemake`
   * `which dot`

