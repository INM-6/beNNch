# Continuous benchmarking

## Benchmark configuration

### Benchmark definitinon
```yaml
repositories:
  [REPOSITORY_ID]:
    url: ""  # 
    branch: ""  # 
benchmarks:
  [MODEL_ID]:
    repository: [REPOSITORY_ID]
    model_dir: ""  # path to model directory inside repository
    machine_parameters:
      [MACHINE_ID]:
        num_nodes: ""  # comma seperated list of integers
        tasks_per_node: ""
        threads_per_task: ""
        walltime: "HH:MM:SS"
    model_parameters:
      - {name: "", type: "float|int|string|text", _: ""}  # 
    run_parameters:
      log_path: ""
      script_path: ""
      run_args: ""
      optional_run_command: ""
    results: []  # list of strings
    metrics: []
    metrics_sum: []
user_config:
  [MACHINE_ID]:
    build-time-limit: "HH:MM:SS"  # slurm time limit for the build job
dependencies: ""  # path to a yaml file with dependency specifications (see #Benchmark dependencies)
```

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

### Model templating


### Automatic result upload


## Pipeline configuration

### Adding a new machine


### Adding a new runner type

