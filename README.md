# Continuous benchmarking

## Benchmark configuration

### Benchmark definitinon
```yaml
repositories:
  common:
    url: ""  # 
    branch: ""  # 
benchmarks:
  hpc_benchmark:
    repository: common
    model_dir: hpc_benchmark  # path to model directory inside repository
    machine_parameters:
      jureca:
        - {name: num_nodes, type: int, _: "1"}
        - {name: tasks_per_node, type: int, _: "2"}
        - {name: threads_per_task, type: int, _: "64"}
        - {name: walltime, type: string, _: "00:02:30"}
        - {name: affinity, type: string, separator: ";", _: "--cpu-bind=verbose,threads --distribution=block:cyclic:fcyclic --threads-per-core=1"}  # processor affinity/pinning
      jusuf:
        - {name: num_nodes, type: int, _: ""}
        - {name: tasks_per_node, type: int, _: ""}
        - {name: threads_per_task, type: int, _: ""}
        - {name: walltime, type: string, _: "00:02:30"}
        - {name: affinity, type: string, separator: ";", _: "--cpu-bind=verbose,threads --distribution=block:cyclic:fcyclic --threads-per-core=1"}  # processor affinity/pinning
    model_parameters:
      - {name: "", type: "float|int|string|text", _: ""}  # 
    run_parameters:
      script_path: ""
    results: []
    metrics: []
    metrics_sum: []
user_config:
  jureca:
    build-time-limit: "00:10:00"  # slurm time limit for the build job
  jusuf:
    build-time-limit: "00:10:00"
dependencies: dependencies.yml
```

### Benchmark dependencies


### Model templating


### Automatic result upload


## Pipeline configuration

### Adding a new machine


### Adding a new runner type

