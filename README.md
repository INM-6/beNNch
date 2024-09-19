# Continuous benchmarking

## Benchmark configuration

### Benchmark definitinon
```yaml
repositories:
  [REPOSITORY_ID]:  # Replace the REPOSITORY_ID by a custom unique identifier for the repository.
    url: ""  # URL to a git repository in a form that makes it possible to run 'git clone [url]'.
    branch: ""  # Branch name of any branch in the repository.
software_parameters:
  modulefile: ""  # Path to a modulefile with instructions to load the simulator module once built
  planfile: ""  # Path to a builder planfile with instructions to get and build the simulator
  simulator_repo_url: ""  # Git repository url which will be forwarded to the builder planfile
  simulator_repo_branch: ""  # Git repository branch which will be forwarded to the builder planfile
benchmarks:
  [MODEL_ID]:
    repository: [REPOSITORY_ID]
    model_dir: ""  # Path to model directory inside repository.
    machine_parameters:
      [MACHINE_ID]:  # Replace the MACHINE_ID by a custom unique identifier for the machine. Must match the machine identifier used to define the machine-specific subpipeline (see #Pipeline configuration).
        num_nodes: ""  # Passed to slurm to set the number of compute nodes to use for the benchmark job.
        tasks_per_node: ""  # Passed to slurm to set the number of MPI processes per compute node to use for the benchmark job.
        threads_per_task: ""  # Passed to slurm to set the number of threads per MPI process to use for the benchmark job.
        walltime: "HH:MM:SS"  # Passed to slurm to set the upper time limit for the benchmark job.
    model_parameters:
      - {name: "", type: "float|int|string|text", _: ""}  # Any parameter that should be substituted in the model 'script_path' file (see #Model templating).
    run_parameters:
      run_command: ""  # Which command to execute to run the script, e.g. 'python'.
      script_path: ""  # Path to the Python script which will be run to start the benchmark for the specific model.
      run_args: ""  # Any arguments that should be passed to the script. Can contain variables using $VARIABLE which come either from the "model_parameters" or the globally defined variables (see # Global variables).
      log_path: ""  # Path where the model stores logfiles which are used to extract timer data.
      optional_run_command: ""  # Any command or list of commands that should be run prior to the script given via the "script_path".
    params_rec: []  # List of strings containing all model parameter keys which should be added as a column to the result.csv file
    metrics_mean: []  # List of strings containing all output parameters in log files which should be added as a column to the result.csv file. Aggregation function for multiple MPI processes (multiple log files): Average of all values
    metrics_sum: []  # List of strings containing all output parameters in log files which should be added as a column to the result.csv file. Aggregation function for multiple MPI processes (multiple log files): Sum of all values
    metrics_step: []  # List of strings containing all output parameters that should be recorded every simulation step.
user_config:
  [MACHINE_ID]:  # Replace the MACHINE_ID by a custom unique identifier for the machine. Must match the machine identifier used to define the machine-specific subpipeline (see #Pipeline configuration).
    build-time-limit: "HH:MM:SS"  # Upper time limit for the build job passed to slurm.
dependencies: ""  # Path to a yaml file with dependency specifications (see #Benchmark dependencies).
```
TODO: Provide yaml-schema for validation.

Any parameter in the *machine_parameters*, *model_parameters*, or *run_parameters* section can receive multiple values per key (e.g. `num_nodes: "1,2,4,8"`). One benchmark run will be executed per possible parameter combination.  

The *model_parameters* and *machine_parameters* can also depend on other parameters or include any arbitrary python code to generate the value or list of values. In place of the value or list of values, simply write python code returning the desired values. For the *model_parameters*, one must additionally add `mode: python`. E.g.:
`{name: scale, mode: python, type: float, "_": "$scale_N*$num_nodes"}`
or
`threads_per_task: "128/$tasks_per_node"`.

Examples can be found here:
https://github.com/JanVogelsang/CI-Bennch-config/blob/deliver-events-first/config.yml

Variables are evaluated right when the pipeline is created (i.e., when clicking on the "Run pipeline" button). If a variable is only known on a specific machine and should be evaluated once jobs of the pipeline run on the machine, use double dollar signs (`$$`) instead of a single one inside the _variables_ section to escape the variable (e.g. `TEMP_DIR: "$$SCRATCH_cjinb33"`). Additionally, any variable containing an escaped variable must be evaluated at runtime (e.g. `'eval TEMP_DIR="$TEMP_DIR"'`) in both the _init_ and _load_ sections.

#### Benchmark dependencies
```yaml
[MACHINE_ID]:
  variables: []  # Any variables that should be available in "vars", "init", and "load"
  vars: []  # Preparing the environment variables. Called before "init" and "load"
  init: []  # Called once to initialize the environment. Called after "vars" but before "load".
  load: []  # Called each time after calling "vars" (but never before "init" was called). Used to load the environment.
```

### Pipeline variables
There is a set of variables beNNch provides which can be used at various places. Some are set by beNNch and some must be
provided by the user. All of the variables below can be safely accessed by the user at any time, given the variable is
already set at that time. The columns 1 and 2 indicate two different points in time when the user might access
these variables. An "X" indicates the variable has already been set at this point.  
1: In the "run pipeline" form
2: In simulator-modulefile/-planfile or dependencies

| Variable name       | Type                            | Read-only | Mandatory | 1 | 2 | Description                                                                                                                                                                                            |
|---------------------|---------------------------------|-----------|-----------|---|---|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| BENCHMARK_ID        | String                          | No        | Yes       | X | X | Label given to the CI pipeline to identify a benchmark.                                                                                                                                                |
| BENCHMARK_UID       | String (POSIX compatible)       | Yes       | /         |   | X | Label converted to a POSIX-compatible version for filenames.                                                                                                                                           |
| CONFIG_REPO         | Any valid input to git clone    | No        | Yes       | X | X | Configuration repository to use for benchmarking.                                                                                                                                                      |
| CONFIG_BRANCH       | String                          | No        | Yes       | X | X | Configuration repository branch to use for benchmarking.                                                                                                                                               |
| CONFIG_PATH         | Relative path to file           | No        | Yes       | X | X | Path to configuration file (e.g. config.yml) relative to the config repository root.                                                                                                                   |
| RESULTS_REPO_DIR    | Absolute path to directory      | No        | Yes       | X | X | Directory of git-annex-enabled repository where result files are moved to.                                                                                                                             |
| EXTRA_CMAKE_FLAGS   | Any valid input to cmake        | No        | No        | X | X | Any extra CMake flags which should be used to compile NEST. Only single quotes allowed! Multiple values for a single parameter must be single quoted and separated by semicolons.                      |
| CB_MODELS           | Comma-separated list of strings | No        | No        | X | X | Comma-separated list of the names of all models for which benchmarks should be run. Leave empty to run benchmarks for all models specified in the configuration file.                                  |
| CB_MACHINES         | Comma-separated list of strings | No        | No        | X | X | Comma-separated list of the names of all machines on which benchmarks should be run. Leave empty to trigger benchmark runs for manually for each individual machine for which configuration was found. |
| CB_SIMULATOR_REPO   | URL to git repository           | Yes       | /         |   | X | The url to the simulator git repository, which will most likely be used as input to git clone in the builder planfile.                                                                                 |
| CB_SIMULATOR_BRANCH | String                          | Yes       | /         |   | X | The simulator git branch to be used.                                                                                                                                                                   |

### Model templating
Find examples for MAM, Microcircuit and HPC benchmark here: 
https://github.com/JanVogelsang/ci-bennch-models.git

### Automatic result upload
In the future, results will be uploaded to git-annex automatically including plenty of metadata (see branch ci-bennch-dev).  
Currently, the data has to be extracted manually. By default, the output data is located under `/p/scratch/cjinb33/ci-bennch/jusuf/userdata/$USER/$BENCHMARK_ID/benchmark_results/`.

For debugging your runs, you might want to look at the log files of each job.  
All output any of the jobs produce is located unter `/p/scratch/cjinb33/ci-bennch/jusuf/userdata/$USER/$BENCHMARK_ID/jube_results/` with the following directory structure:
`jube_results`: This folder is completely managed by Jube and contains one subdirectory per Jube-Run.
`jube_results/000000` (e.g.): In ci-bennch, one Jube-Run is started per model that is being benchmarked. For each model, one slurm job is started for each parameter combination.
`jube_results/000000/000000_bench` (e.g.): This directory contains all files required for a specific job with all template parameters replaced by actual values for a specific parameter combination. It also contains any output data, logs, and files created during the job.

## Pipeline configuration

### Adding a new machine


### Adding a new runner type

