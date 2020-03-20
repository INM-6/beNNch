
## The control environment

   In the descriptions of the benchmarking workflow, the most outside
   environment you are in to start the benchmarks will be called the *control
   environment*.  This is usually the shell from which you start everything on
   a machine sufficiently close to where you want to run the benchmarks (e.g.
   the head node of a cluster).

   The control environment basically has to provide only conda to bootstrap the
   benchmark workflow. You may, however, have a version of miniconda installed
   already (ref.  [install]()) with snakemake in your base environment, and
   most machines provide `wget` by default. In this case you are all set and do
   not explicitly need to set up an extra control environment.

   The following sections detail some choices you have to make your life
   easier.


### Shell alias

   All calls of the workflow are made via snakemake and a couple of options
   should be given every time. It it handy to create a shell alias to not type
   all of them each time and potentially forget one every now and then.

   To make your life easier add the following line to your shell environment
   (e.g.  in `~/.bashrc`):

       alias smake='snakemake -j --use-conda'

   Having this available (after reopening the shell) you can just run `smake`
   instead of `snakemake`.


### Conda environment locations

   Conda environments can be voluminous you may want to keep them in different
   folder than the default ".snakemake/conda" inside your workflow repository
   (which is where Snakemake puts them by default).

   To change this you can give `snakemake` a directory setting `--conda-prefix`
   (see [Integrated Package Management][] in the Snakemake docs) to define
   where to store the automatically created environments.

    snakemake --use-conda --conda-prefix some/path/to/storage

   It is easiest to include this flag in your [Shell alias](#shell-alias).

   [Integrated Package Management]: https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management Snakemake documentation – Integrated Package Management


### Explicit control environment

   In case you are on a very restricted machine, your control environment, may
   well be a conda environment set up with the `environment/control.yaml`, which
   provides all required tools.

       conda env create --name bench --file environment/control.yaml
       conda activate bench

   To run the benchmark workflow be sure to be in this environment first, then
   no additional setup is required. Note however, that this may not be
   necessary if tools can be provided by other means, e.g. by
   [modules](#about-modules).


### About modules

   In cluster environments software can be made available by *modules*, which
   can be listed with `module avail` command. Some tools like `wget` may for
   example be available only after running

       module load wget

   Using the available tools may lift the requirement to set up an [explicit
   control environment](#explicit-control-environment)

