# CI-beNNch

[![License: GPL v3+](https://img.shields.io/badge/License-GPLv3%2B-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Cite this repository](https://img.shields.io/badge/cite-CITATION.cff-green.svg)](CITATION.cff)

**A modular framework for performance benchmarking of neuronal network simulations.** CI-beNNch makes benchmark runs of simulation engines such as [NEST](https://www.nest-simulator.org/) reproducible and comparable across machines for both early-career and highly experienced researchers.

## Status and scope

CI-beNNch has been in operation for roughly four years at Forschungszentrum Jülich, on Jülich Supercomputing Centre (JSC) supercomputers and on institute-local clusters. This repository is the **public-facing** home of CI-beNNch. The pipeline definitions target GitLab CI with Jacamar CI runners attached to JSC compute resources. Those runners live in a private GitLab instance (`gitlab.jsc.fz-juelich.de`) that is not reachable without a JSC account, and the benchmark jobs require a JSC compute-time allocation. You therefore cannot run this repository as-is from GitHub. To run continuous benchmarking on your own infrastructure, use this code as a reference implementation. The structure of the pipeline, the machine configurations and the generation scripts are the reusable parts.

## Quickstart

Running CI-beNNch means pointing a GitLab pipeline at a benchmark configuration.

**Prerequisites:** a GitLab instance with Jacamar CI runners registered against your HPC system with a compute budget on that system.

1. **Mirror this repository** to your GitLab instance. `.gitlab-ci.yml` at the root is the pipeline entry point.

2. **Describe your machine.** Copy the closest match from [`.gitlab-ci/machine-configs/`](.gitlab-ci/machine-configs/) (for example [`jureca.yml`](.gitlab-ci/machine-configs/jureca.yml)) and adjust the runner tags, partitions and module setup. The file name (minus `.yml`) becomes the `MACHINE_ID` you use everywhere else.

3. **Write a benchmark configuration** in a separate git repository, see [Benchmark configuration](#benchmark-configuration) below for the full schema. Keeping it separate from the pipeline code is deliberate.

4. **Run the pipeline.** In GitLab, go to *Build → Pipelines → Run pipeline*, and set at least `BENCHMARK_ID`, `CONFIG_REPO`, `CONFIG_BRANCH` and `CONFIG_PATH` (see [Pipeline variables](#pipeline-variables)). Pipelines only start when triggered from the web UI, by a trigger token, or on a schedule. Pushes do not automatically start benchmarks.

5. **Collect results.** Output lands under the path described in [Result handling](#result-handling).

## How to cite

If you use CI-beNNch in published work, please cite it. GitHub's *"Cite this repository"* button reads [`CITATION.cff`](CITATION.cff) and will give you a ready-made entry.

For version 2.0 (CI-beNNch):

> Vogelsang, J., Lober, M., Schöfmann, C. M., Villamar, J., Terhorst, D., Senk, J.,
> Plesser, H. E., Diesmann, M., Kunkel, S., & Kurth, A. C. (2026). Continuous
> benchmarking: Keeping pace with an evolving ecosystem of models and
> technologies. *arXiv preprint* arXiv:2604.15919.
> doi: [10.48550/arXiv.2604.15919](https://doi.org/10.48550/arXiv.2604.15919)

For the original workflow (beNNch 1.x):

> Albers, J., Pronold, J., Kurth, A. C., Vennemo, S. B., Haghighi Mood, K.,
> Patronis, A., Terhorst, D., Jordan, J., Kunkel, S., Tetzlaff, T., Diesmann, M.,
> & Senk, J. (2022). A Modular Workflow for Performance Benchmarking of Neuronal
> Network Simulations. *Frontiers in Neuroinformatics*, 16, 837549.
> doi: [10.3389/fninf.2022.837549](https://doi.org/10.3389/fninf.2022.837549)

Machine-readable metadata: [`CITATION.cff`](CITATION.cff) and [`codemeta.json`](codemeta.json).

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for how to report a bug, what to include, and how changes are reviewed.

## License

GNU General Public License v3.0 or later (`GPL-3.0-or-later`) (see [LICENSE](LICENSE)). Copyright © Forschungszentrum Jülich GmbH, IAS-6.

---

# Reference

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

*Not yet implemented: a YAML schema for validating this configuration.*

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

## Pipeline variables

### Variables you set when starting a pipeline

Defined in [`.gitlab-ci.yml`](.gitlab-ci.yml); these appear as fields in the GitLab *Run pipeline* form.

| Variable name | Mandatory | Description |
|---|---|---|
| `BENCHMARK_ID` | Yes | Label given to the CI pipeline to identify a benchmark. Becomes the pipeline name. |
| `CONFIG_REPO` | Yes | Configuration repository to use for benchmarking (anything `git clone` accepts). |
| `CONFIG_BRANCH` | Yes | Configuration repository branch to use. |
| `CONFIG_PATH` | Yes | Path to the configuration file (e.g. `config.yml`) relative to the config repository root. |
| `RESULTS_REPO_URL` | Yes | Repository to upload benchmark results to. |
| `RESULTS_BRANCH` | Yes | Branch of the result repository to upload results to. |
| `RESULTS_PATH` | Yes | Path to the results directory relative to the results repository root. |
| `CB_MODELS` | No | Model(s) for which benchmarks should be run, e.g. `microcircuit`, `hpc_benchmark`, `multi_area`, `cmam`, `cmam_connect_benchmark`, `mesocircuit`. Leave empty to run all models in the configuration file. |
| `CB_MACHINES` | No | Machine(s) on which benchmarks should be run, e.g. `jureca`, `juwels`, `jusuf`, `jupiter`, `hambach`. Leave empty to trigger each machine's run manually. |
| `CB_STORAGE_MACHINE` | No | Runner tags of the storage machine. |
| `CB_STORAGE_DIR` | No | Directory of the repository on `CB_STORAGE_MACHINE`. |
| `EXTRA_CMAKE_FLAGS` | No | Extra CMake flags for compiling NEST. Single quotes only; multiple values for one parameter must be single-quoted and semicolon-separated. |
| `EXTRA_SBATCH_FLAGS` | No | Extra `sbatch` flags used when running on the machine. |
| `DETAILED_TIMERS_FLAG` | No | `ON`/`OFF` - whether detailed timers are shown. |
| `MPI_FLAG` | No | `ON`/`OFF` - whether NEST is compiled with MPI. |
| `PLOTTING_STYLE` | No | `Standard`/`Docs`. `Docs` forces detailed timers `OFF` and plots the realtime factor instead. |

### Variables beNNch derives

Set by beNNch itself; read-only. Available in the dependencies file.

| Variable name | Description |
|---|---|
| `BENCHMARK_UID` | `BENCHMARK_ID` converted to a POSIX-compatible form, for use in filenames. |
| `CB_SIMULATOR_REPO` | URL of the simulator git repository. |
| `CB_SIMULATOR_BRANCH` | Simulator git branch being used. |

## Model templating

Example model templates for the multi-area model, the microcircuit and the HPC benchmark are maintained alongside this project. See the release notes for the location current at the time of release.

## Result handling

By default the output data is written to `/p/scratch/cjinb33/ci-bennch/<machine>/userdata/$USER/$BENCHMARK_ID/benchmark_results/` on the JSC systems.

For debugging runs, the log files of each job are located under `/p/scratch/cjinb33/ci-bennch/<machine>/userdata/$USER/$BENCHMARK_ID/jube_results/` with the following structure:

- `jube_results` - managed entirely by JUBE; contains one subdirectory per JUBE run.
- `jube_results/000000` - one JUBE run is started per benchmarked model. For each
  model, one Slurm job is started per parameter combination.
- `jube_results/000000/000000_bench` - all files required for a specific job, with
  template parameters replaced by the actual values for that parameter
  combination, plus any output data, logs and files created during the job.
