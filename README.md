# NEST benchmarks

### Current repository layout

*analysis* contains JUBE analysis script, config and helpers.

*benchmarks* contains benchmark scripts to run benchmarks via JUBE.

*config* contains user configuration file templates to be copied and adapted. 

*helpers* contains helper JUBE parametersets.

*results* contains all results and analysis scripts.


### Initialization

- download submodule data
  + `git submodule init`
  + `git submodule update`
- install benchplot as module
  + `pip install -e plot --user`

### First steps

Copy `config/templates/user_config_template.xml` to `config/user_config.xml` and fill in all parameters:
  - `account`: slurm account name for job submission
  - `email_address` (optional): email address to which slurm sends START, END, FAIL emails for job progress info
  - `partition` (optional, required on some machines): cluster partition, takes system default if left empty 

Copy `config/templates/<model>_config_template.xml` to `config/<model>_config.xml` and fill in all parameters.


### Software Dependencies

- [git annex](https://git-annex.branchable.com)
  + can be installed via 
```bash
wget 'http://downloads.kitenet.net/git-annex/linux/current/git-annex-standalone-amd64.tar.gz'
tar -xzf git-annex-standalone-amd64.tar.gz
export PATH=$PATH:<install_path>/git-annex.linux
```
- [JUBE](https://www.fz-juelich.de/ias/jsc/EN/Expertise/Support/Software/JUBE/_node.html)
- Python 3


### Running Benchmarks

The JUBE scripts for the benchmarks can be found in `benchmarks/`.

To run a benchmark, run

```bash
jube run <PATH_TO_REPO>/benchmarks/<benchmark_file.xml>
```

You will get a job `id`.

These are the benchmarks currently implemented:

- **Multi-Area Model**

  - `benchmarks/multi-area-model_2.xml` for usage with NEST 2.
  - `benchmarks/multi-area-model_3.xml` for usage with NEST 3.

- **Microcircuit**

  - `benchmarks/microcircuit.xml`

- **HPC Benchmark**

  - `benchmarks/hpc_benchmark_2.xml` for usage with NEST 2.
  - `benchmarks/hpc_benchmark_3.xml` for usage with NEST 3.0.
  - `benchmarks/hpc_benchmark_31.xml` for usage with NEST 3.1.

### Analyze Benchmarks

- copy `analysis/analysis_config_template.py` to `analysis/analysis_config.py`
- fill in
  + type of scaling (for creating a quick, glanceable plot of the benchmark. Here we provide defaults for plotting timers across `nodes` and `threads`. To create your own plot, add to `analysis/plot_helpers.py`.)
  + path to the jube output (usually the same as the `outpath` of the `<benchmark>` in `benchmarks/<model>`)
- `cd results` (s.t. git annex metadata annotation works)
- `python ../analysis/analysis.py <id>` where `<id>` is the JUBE ID of the benchmark you want to analyze
  + if this is the first time `results` is used, get up to date via `git pull origin main`, `git checkout main`, `git annex sync`
- if you're happy with the results: `git annex sync`

### Get remote benchmark results
- `cd results`
- `git annex init` (only needed at first usage for initialisation)
- `git add remote <name> <location>`, e.g. `git add remote jusuf /p/project/icei-hbp-2020-0006/ACA_bm_framework/nest-benchmarks/results`
- `git fetch <name>`
- `git annex get`

### Visualization

- go to `results`
- select benchmarks to plot via the following syntax:
  + `git annex view <common_metadata>=<value_of_metadata> <differing_metadata>="*"`
    * Here, common_metadata refers to a key that should be the same value for all benchmarks, e.g. the `machine`. A list of all available metadata keys can be obtained via `git annex metadata <uuidgen_hash>.csv`. One can use `*` here as well, e.g. when filtering out all runs that include simulations done on 10 nodes via `num_nodes='*,10*'` or all machines that have `jusuf` in their name via `machine='*jusuf*`. To specify multiple numerical values use `keyword={value1,value2}`.
    * full example: `git annex view nest=nest-simulator/3.0 num_vps="*"`
    * to go back in a view, execute `git annex vpop`
- create slideshow of plots with `python ../slideshow/slideshow.py <scaling_type> <bullet_1> <bullet_2> ...` with an arbitrarily long list of bullet items (metadata keys) that appear as bullet points on the slides for comparison. `<scaling_type>` defines the style of plotting, c.f. section on _Analyze Benchmarks_.

#### Known issues
- error `jinja2.exceptions.TemplateNotFound: index.html.j2`
  + [issue](https://github.com/jupyter/nbconvert/issues/1394) with a recent version of `nbconvert`, try to install version `5.6.1` instead (e.g. `pip install nbconvert==5.6.1 --user`)


#### Some benchmark information

**Population model**

The model has several populations (minimum 20 populations). Each population has 5000 neurons, and each neuron has 5000 connections. Each population connects to 20 other populations, which are chosen randomly. There are thus 250 connections per neuron and projection. The way we scale up is by adding more populations. We still connect to 20 other populations.

#### To run MAM:

You need to have installed NEST *with* Python in order to run the MAM benchmark.

You also might need to download *nested_dict* and *dicthash*:

```bash
pip install nested_dict --user
pip install dicthash --user
```
















