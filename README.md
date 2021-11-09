# NEST Benchmarking Framework

## Repository structure

*analysis* contains JUBE analysis script, config and helpers.

*benchmarks* contains benchmark scripts to run benchmarks via JUBE.

*config* contains user configuration file templates to be copied and adapted. 

*helpers* contains helper JUBE parametersets.

*results* contains all results and analysis scripts.

*models* is a git submodule; the linked repository contains NEST network models adapted to work with the benchmarking framework.

*plot* is a git submodule; the linked repository contains predefined plotting routines designed to process the performance results and provide a standardized plotting format.


## Using the framework

### Initialization

- download submodule data
  + `git submodule init`
  + `git submodule update`
- install benchplot as module
  + `pip install -e plot --user`
  
### Software dependencies

- [git annex](https://git-annex.branchable.com)
  + can e.g. be installed via 
```bash
wget 'http://downloads.kitenet.net/git-annex/linux/current/git-annex-standalone-amd64.tar.gz'
tar -xzf git-annex-standalone-amd64.tar.gz
export PATH=$PATH:<install_path>/git-annex.linux
```
- [JUBE](https://www.fz-juelich.de/ias/jsc/EN/Expertise/Support/Software/JUBE/_node.html)
- Python 3.X

### First steps: configure your simulation

Copy `config/templates/user_config_template.yaml` to `config/user_config.yaml` and fill in all parameters:
  - `account`: slurm account name for job submission
  - `email_address` (optional): email address to which slurm sends START, END, FAIL emails for job progress info
  - `partition` (optional, required on some machines): cluster partition, takes system default if left empty 

Copy `config/templates/<model>_config_template.yaml` to `config/<model>_config.yaml` and fill in all parameters.

### Run benchmarks

The JUBE benchmarking scripts can be found in `benchmarks/`.

To run a benchmark, execute

```bash
export JUBE_INCLUDE_PATH="config/:helpers/"
jube run <PATH_TO_REPO>/benchmarks/<benchmark_file.yaml>
```

_Note that the export command is working around a known bug in JUBE 2.4.1, once it is fixed the export will become unnecessary and the documentation here will be updated accordingly._

JUBE displays a table summarizing the submitted job(s) and the corresponding `job id`.

These are the benchmarks currently implemented:

- **Multi-Area Model**

  - `benchmarks/multi-area-model_2.yaml` for usage with NEST 2.
  - `benchmarks/multi-area-model_3.yaml` for usage with NEST 3.

- **Microcircuit**

  - `benchmarks/microcircuit.yaml`

- **HPC Benchmark**

  - `benchmarks/hpc_benchmark_2.yaml` for usage with NEST 2.
  - `benchmarks/hpc_benchmark_3.yaml` for usage with NEST 3.0.
  - `benchmarks/hpc_benchmark_31.yaml` for usage with NEST 3.1.

### Analyze benchmarks

First, create a new instance of the analysis configuration with
```bash
cp analysis/analysis_config_template.py analysis/analysis_config.py
```
Here, fill in
- whether the scaling benchmark runs across threads or nodes. This sets up a quick, glanceable plot of the benchmark to confirm that no substatial errors occurred. The framework provides defaults for plotting timers across `nodes` and `threads`, but alternatives can be readily implemented by adding to `analysis/plot_helpers.py`.
- the path to the JUBE output (usually the same as the `outpath` of the `<benchmark>` in `benchmarks/<model>`)

To start the analysis, execute
```bash
cd results
```
- _optional: initialize for the first time_
  + `git pull origin main`
  + `git checkout main`
  + `git annex init`
  + `git annex sync`
```bash
python ../analysis/analysis.py <id>
```
where `<id>` is the `job id` of the benchmark you want to analyze.

For sharing, upload the results to the central repository via
```bash
git annex sync
```

### Get remote benchmark results
```bash
cd results
```
- _optional: add a new remote_ 
  + `git add remote <name> <location>`, e.g. `git add remote jureca <username>@jureca.fz-juelich.de:<PATH/TO/REPO>/results`
  + `git fetch <name>`
```bash
git annex get
```

### Visualization

```bash
cd results
```

First, filter which benchmarks to plot using the following syntax:
```bash
git annex view <common_metadata>="<value_of_metadata>" <differing_metadata>="*"
```
Here, `common_metadata` refers to a key that should be the same value for all benchmarks, e.g. the `"machine"`. A list of all available metadata keys can be obtained via `git annex metadata <uuidgen_hash>.csv`. One can use `*` here as well, e.g. when filtering out all runs that include simulations done on 10 nodes via `num_nodes='*,10*'` or all machines that have `jureca` in their name via `machine='*jureca*`. To specify multiple numerical values use `keyword={value1,value2}`.
- example: `git annex view machine="jureca" model_name="microcircuit" nest="*" num_vps="*"`.

Note that this changes the local file structure; the values corresponding to the `<differing_metadata>` determine the names of the folders in a hierarchical fashion. In the example above, the top level would consist of folders named after the values of `nest` (e.g. `nest-simulator/2.14.1`, `nest-simulator/2.20.2` and `nest-simulator/3.1`), with each of those containing folders named after the number of virtual processes of the simulations (e.g. `4`, `8`, `16`). Rearranging the order of `<differing_metadata>` in the command above also reorders the hierarchical file structure.
To "go back" a view, execute
```bash
git annex vpop
```
After choosing which benchmarks to display via filtering above and ordering them via `<differing_metadata>`, you can create a slideshow of all plots with
```bash
python ../slideshow/slideshow.py <scaling_type> <bullet_1> <bullet_2> ...
```
with an arbitrarily long list of bullet items (consisting of metadata keys) that appear as bullet points on the slides for comparison. `<scaling_type>` defines the style of plotting, c.f. section on [Analyze Benchmarks](#analyze-benchmarks).

### Known issues
- error `jinja2.exceptions.TemplateNotFound: index.html.j2`
  + [issue](https://github.com/jupyter/nbconvert/issues/1394) with a recent version of `nbconvert`, try to install version `5.6.1` instead (e.g. `pip install nbconvert==5.6.1 --user`)

___
## Developing the framework

### Add a new model

`benchmarks/template.yaml` provides a template for a JUBE benchmarking script and can be used as a starting point for adding a new model. Here, only the marked section needs to be adapted. As a reference, see the implementation of the microcircuit in `benchmarks/microcircuit.yaml`.

In addition, minor modifications to a regular network model need to be made in order to comply with the framework's standards. In particular, this concerns how JUBE feeds the configuration parameters to the network and how JUBE reads the performance measurement output.

#### Input

A new model needs to be able to receive input from JUBE for setting parameters. Following the `substituteset` defined in `benchmarks/template.yaml`, all listed source keys need to be initialized. See `models/Potjans_2014/run_bm_microcircuit.py` for reference.

#### Output

As current releases of NEST (including 2.14.1, 2.20.2 and 3.0+) include timers on the C++ level for measuring the simulation performance, the model only needs to output this information in a way compliant with the framework. This can be done via adding a call to the `logging` function defined in `models/Potjans_2014/bm_helpers.py`. Note that this also provides the optional functionality to include python level timers as well as memory information.
