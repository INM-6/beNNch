# CI-beNNch — continuous benchmarking of neuronal network simulations

> This repository is the **public-facing** home of CI-beNNch. It was developed to execute benchmarks as CI pipelines on HPC systems at the Jülich Supercomputing Centre. The pipeline definitions in this repository target GitLab CI with [Jacamar CI](https://gitlab.com/ecp-ci/jacamar-ci) runners attached to JSC compute resources. Those runners live in a private GitLab instance (`gitlab.jsc.fz-juelich.de`) that is not publicly reachable without an JSC account, and the benchmark jobs require a JSC compute-time allocation.
> If you want to run continuous benchmarking on **your own** infrastructure, our code serves as a reference.

## How to cite

For the original workflow and the concepts both versions build on, please cite:

> Albers J, Pronold J, Kurth AC, Vennemo SB, Haghighi Mood K, Patronis A, Terhorst D,
> Jordan J, Kunkel S, Tetzlaff T, Diesmann M and Senk J (2022) A Modular Workflow for
> Performance Benchmarking of Neuronal Network Simulations. *Front. Neuroinform.* 16:837549.
> doi: [10.3389/fninf.2022.837549](https://doi.org/10.3389/fninf.2022.837549)

To cite a specific version, use the release DOI from the Zenodo record for that release, or
name the release and this repository.

## License

GNU General Public License v3.0 — see [LICENSE](LICENSE).

---

# Documentation

## Benchmark configuration

### Benchmark definition

```yaml
repositories:
  <REPOSITORY_ID>:  # Replace the REPOSITORY_ID by a custom unique identifier for the repository.
    url: ""  # URL to a git repository in a form that makes it possible to run 'git clone <url>'.
    branch: ""  # Branch name of any branch in the repository.
software_parameters:
  simulator_repo_url: ""  # Git repository url to clone the simulator from
  simulator_repo_branch: ""  # Git repository branch to use when cloning the simulator
  build_instructions: []  # List of commands to build the simulator
  install_instructions: []  # List of commands to install/load the simulator
benchmarks:
  <MODEL_ID>:
    repository: <REPOSITORY_ID>
    model_dir: ""  # Path to model directory inside repository.
    machine_parameters:
      <MACHINE_ID>:  # Replace the MACHINE_ID by a custom unique identifier for the machine. Must match the machine identifier used to define the machine-specific subpipeline (see #Pipeline configuration).
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
  <MACHINE_ID>:  # Replace the MACHINE_ID by a custom unique identifier for the machine. Must match the machine identifier used to define the machine-specific subpipeline (see #Pipeline configuration).
    build-time-limit: "HH:MM:SS"  # Upper time limit for the build job passed to slurm.
dependencies: ""  # Path to a yaml file with dependency specifications (see #Benchmark dependencies).
```

*Not completed before development ceased: a YAML schema for validating this configuration.*

Any parameter in the *machine_parameters*, *model_parameters*, or *run_parameters* section can receive multiple values per key (e.g. `num_nodes: "1,2,4,8"`). One benchmark run will be executed per possible parameter combination.

The *model_parameters* and *machine_parameters* can also depend on other parameters or include any arbitrary python code to generate the value or list of values. In place of the value or list of values, simply write python code returning the desired values. For the *model_parameters*, one must additionally add `mode: python`. E.g.:
`{name: scale, mode: python, type: float, "_": "$scale_N*$num_nodes"}`
or
`threads_per_task: "128/$tasks_per_node"`.

Variables are evaluated right when the pipeline is created (i.e., when clicking on the "Run pipeline" button). If a variable is only known on a specific machine and should be evaluated once jobs of the pipeline run on the machine, use double dollar signs (`$$`) instead of a single one inside the _variables_ section to escape the variable.

#### Benchmark dependencies

```yaml
<MACHINE_ID>:
  variables: []
```

### Pipeline variables

There is a set of variables beNNch provides which can be used at various places. Some are set by beNNch and some must be
provided by the user. All of the variables below can be safely accessed by the user at any time, given the variable is
already set at that time. The columns 1 and 2 indicate two different points in time when the user might access
these variables. An "X" indicates the variable has already been set at this point.
1: In the "run pipeline" form
2: In the dependencies file

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
| CB_SIMULATOR_REPO   | URL to git repository           | Yes       | /         |   | X | The url to the simulator git repository.                                                                                       |
| CB_SIMULATOR_BRANCH | String                          | Yes       | /         |   | X | The simulator git branch to be used.                                                                                                                                                                   |

### Model templating

Example model templates for the multi-area model, the microcircuit and the HPC benchmark are
maintained alongside this project. See the release notes for the location current at the time
of release.

### Result handling

Automatic upload of results to git-annex with accompanying metadata was designed but **not
completed** before development ceased; in 2.0 the data has to be extracted manually. By
default the output data is written to
`/p/scratch/cjinb33/ci-bennch/<machine>/userdata/$USER/$BENCHMARK_ID/benchmark_results/`
on the JSC systems.

For debugging runs, the log files of each job are located under
`/p/scratch/cjinb33/ci-bennch/<machine>/userdata/$USER/$BENCHMARK_ID/jube_results/` with the
following structure:

- `jube_results` — managed entirely by JUBE; contains one subdirectory per JUBE run.
- `jube_results/000000` — one JUBE run is started per benchmarked model. For each model, one
  Slurm job is started per parameter combination.
- `jube_results/000000/000000_bench` — all files required for a specific job, with template
  parameters replaced by the actual values for that parameter combination, plus any output
  data, logs and files created during the job.
