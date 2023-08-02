# Continuous benchmarking

## Benchmark configuration

### Benchmark definitinon
```yaml
repositories:
  [REPOSITORY_ID]:
    url: ""  # 
    branch: ""  # 
benchmarks:
  hpc_benchmark:
    repository: ""
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
      script_path: ""
    results: []  # list of strings
    metrics: []
    metrics_sum: []
user_config:
  [MACHINE_ID]:
    build-time-limit: "00:10:00"  # slurm time limit for the build job
dependencies: "/path/to/dependencies/file.yml"
```

### Benchmark dependencies


### Model templating


### Automatic result upload


## Pipeline configuration

### Adding a new machine


### Adding a new runner type

