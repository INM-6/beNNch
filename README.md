# Continuous benchmarking

## Benchmark configuration

### Benchmark definitinon
```yaml
repositories:
  [REPOSITORY_ID]:  # Replace the REPOSITORY_ID by a custom unique identifier for the repository.
    url: ""  # URL to a git repository in a form that makes it possible to run 'git clone [url]'.
    branch: ""  # Branch name of any branch in the repository.
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
      log_path: ""  # Path where the model stores logfiles which are used to extract timer data.
      script_path: ""  # Path to the Python script which will be run to start the benchmark for the specific model.
      run_args: ""  # Any arguments that should be passed to the script. Can contain variables using $VARIABLE which come either from the "model_parameters" or the globally defined variables (see # Global variables).
      optional_run_command: ""  # Any command or list of commands that should be run prior to the script given via the "script_path".
    results: []  # List of strings 
    metrics: []  # List of strings 
    metrics_sum: []  # List of strings 
user_config:
  [MACHINE_ID]:  # Replace the MACHINE_ID by a custom unique identifier for the machine. Must match the machine identifier used to define the machine-specific subpipeline (see #Pipeline configuration).
    build-time-limit: "HH:MM:SS"  # Upper time limit for the build job passed to slurm.
dependencies: ""  # Path to a yaml file with dependency specifications (see #Benchmark dependencies).
```

Any parameter in the *machine_parameters*, *model_parameters*, or *run_parameters* section can receive multiple values per key (e.g. `num_nodes: "1,2,4,8"`). One benchmark run will be executed per possible parameter combination.

### Benchmark dependencies
```yaml
[MACHINE_ID]:
  variables:
    TEMP_DIR: ""   # Directory to temporary files accessible from both login and compute nodes
    USERSOFTWARE_DIR: ""
  init: []
  load: []  # list of strings
```
Variables are evaluated when the pipeline is created (i.e., when clicking on the "Run pipeline" button). If a variable is only known on a specific machine and should be evaluated once jobs of the pipeline run on the machine, use double dollar signs (`$$`) instead of a single one inside the _variables_ section to escape the variable (e.g. `TEMP_DIR: "$$SCRATCH_cjinb33"`). Additionally, any variable containing an escaped variable must be evaluated at runtime (e.g. `'eval TEMP_DIR="$TEMP_DIR"'`) in both the _init_ and _load_ sections.

### Global variables

### Model templating


### Automatic result upload


## Pipeline configuration

### Adding a new machine


### Adding a new runner type

