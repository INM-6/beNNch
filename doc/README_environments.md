
# Environments

To have the required software tools available for the full chain of
benchmarking scripts, a modular approach has been established. At each level of
detail only a certain subset of tools needs to be available and may depend on
parameters coming from outer layers. The whole framework is built in a way that
each environment creates the environments for the things done as necessary,
abstracting away internal details. As an example, the outermost environment
just needs `conda` and `Snakemake` installed and will create the environment
for the job submission with JUBE on demand. JUBE will in turn create some
environments to build software and eventually other environments to run the
software.

## Control environment

   snakemake is used to run the workflow (local or on head node)

   On this level control of the abstract workflow is organized, defining what
   should be done.

## Launch environment

   JUBE is used to submit jobs to the batch system (usually on the head node)

   Here the interface to the cluster infrastructure is defined. Batch-system
   specific configuration is instanciated and information about the machine and
   job configuration is gathered/generated.

## Job environment

   This is where the submitted job-script runs (usually on a compute node)

   Here the interface to the machine software setup is defined.

   It may in turn set up an environment for the benchmark using `conda` or
   `module load`, instanciating the machine specific environment for jobs.
   From the job environment usually a call is made to `mpirun`, `srun` or one
   of their siblings to start the run.

## Run environment

   Here the actual benchmark runs with everything set up (usually in parallel
   context on all allocated nodes).

   The complete software environment should be available and only model
   specific commands should be necessary here, e.g.  `python
   microcircuit_bench.py`.



# Compile environments

   When compiling the simulator and/or other software, it is important to
   consider how often an environment will be used. For the `launch` environment
   it is very likely only one will be used for all jobs, very few tools just
   need to be there. For the `compile` and `run` environments dependencies and
   installs require a bit more thinking.

   * installs into the environment may collide with other jobs that may run in
     parallel but with different compile time parameters.
     * When using Conda, this workflow uses `conda create --clone` for the
       derived environment (hardlinks, man inodes!). An alternative would be to
       create a new empty environment for the compile and run "conda activate
       --stack" to combine a base and specialized environment.

   * installs not into the environment needs to keep track of this dependency,
     since the compiled software depends on a certain set of libraries in
     specific versions.

   * the environment for running the software needs to provide at least
     *exactly* the same run-time dependencies. It may drop compile-time
     dependencies (like compilers, etc.) and it may add additional dependencies
     (like scipy).

