# NEST benchmarks

To run the benchmarks you need to clone the nest-benchmarks repository. To clone with git, do `git clone https://gin.g-node.org/jasperalbers/nest-benchmarks.git`

This will create the folder `nest-benchmarks`, and it contains all JUBE files needed to run the benchmarks.

### Current repository layout

*benchmarks* contains benchmark scripts to run benchmarks via JUBE.

*helpers* contains helper JUBE parametersets.

*config* contains user configuration file templates to be copied and adapted. 

*results* contains all results and analysis scripts.

### First step

Copy `user_config.xml` to `user_config_<your_ID_name>` and fill in all parameters:
  - `model_path`: path to the neuroscience model
  - `data_path`: path where simulation (spiking) output will be stored
  - `account`: slurm account name for job submission
  - `email_address` (optional): email address to which slurm sends START, END, FAIL emails for job progress info
  - `partition` (optional, required on some machines): cluster partition, takes system default if left empty 

The benchmarks are run using the automatic benchmarking environment [JUBE](https://www.fz-juelich.de/ias/jsc/EN/Expertise/Support/Software/JUBE/_node.html), so if you don't already have this installed, you will need to install it.

### Software Installation

Before running the Benchmarks, you need to install NEST, and you should also install jemalloc.

#### NEST

```bash
cd <PATH>/BenchWork
mkdir NEST
cd NEST
git clone https://github.com/nest/nest-simulator.git src
```

You need to clone into`src`. If you want a specific version of NEST, you need to adjust accordingly, you can find our build files [here](https://gin.g-node.org/nest/nest-benchmarks/src/master/Benchmarks/jube_build) to see which build files we have already created.

Once the source code is in place, build and install NEST via JUBE:

```bash
jube run <PATH>/nest-benchmarks/Benchmarks/jube_build/build_nest_daint_strict.xml
```

NEST will be installed in `BenchWork/NEST-daint-base_O3_strict`.

**NEST with Python**

If you want to build and install NEST with Python, run

```bash
jube run <PATH>/nest-benchmarks/Benchmarks/jube_build/build_nest_py_daint_strict.xml
```

**NEST with Boost**

If you want to run NEST with Boost, you need to manually copy out lines 35, 38, 188,193, 194 and 195 in `<PATH>/BenchWork/NEST/src/libnestutil/sort.h`

![](https://gin.g-node.org/nest/nest-benchmarks/src/master/boost1.png)

![](https://gin.g-node.org/nest/nest-benchmarks/src/master/boost2.png)


#### Jemalloc

```bash
cd <PATH>/BenchWork
mkdir jemalloc
cd jemalloc
wget https://github.com/jemalloc/jemalloc/releases/download/5.0.1/jemalloc-5.0.1.tar.bz2
tar xvf jemalloc-5.0.1.tar.bz2
```

To build and install jemalloc, run

```bash
jube run <PATH>/nest-benchmarks/Benchmarks/jube_build/build_jemalloc_daint.xml
```

### Running Benchmarks

The JUBE scripts for the benchmarks can be found in `<PATH>/nest-benchmarks/benchmarks/`.

To run a benchmark, run

```bash
jube run <PATH>/nest-benchmarks/benchmarks/<benchmark_file.xml>
```

You will get a job `id`.

These are the benchmarks we currently have, with corresponding JUBE files:

- **Multi-Area Model**

  - `multi-area-model_2.xml` for usage with NEST 2.
  - `multi-area-model_2.xml` for usage with NEST 3.

#### legacy models:

- **HPC_benchmark**

  - `hpc_benchmark.xml`

- **HPC_benchmark with static synapse**

  - `hpc_benchmark_daint_strict.xml`

    - where we change the parameter `PLASTIC` to false:

      
        `<parameter name="PLASTIC" type="string">false</parameter>`

- **HPC_benchmark with random delay**

  - `hpc_benchmark_daint_strict.xml`

    - where we change the parameters `D_MIN` and `D_MAX`:

      
        `<parameter name="D_MIN" type="float">0.1</parameter>` 

        `<parameter name="D_MAX" type="float">50.</parameter>`

- **HPC benchmark with different connection rules**

  - `hpc_benchmark_daint_strict_rule.xml`

    - use `in`, `out`, `all`, `one`, `tot`, `bern` or`sym_bern` in `RULE` parameter

      
        `<parameter name="RULE" type="string">in</parameter>`

    - Probably need to use a lower basis_scale to get some of them to work, like `basis_scale=5

- **Population model**

  - `population_daint_strict.xml`

- **Population python model**

  - `population_py_daint_strict.xml`




- **HPC Benchmark on fixed VPs, changing threads**

  - `hpc_benchmark_daint_strict.xml`
    - with `<parameter name="THREADS_PER_TASK" type="int">1,3,6,9,18,36</parameter>`


#### Some benchmark information

**Population model**

The model has several populations (minimum 20 populations). Each population has 5000 neurons, and each neuron has 5000 connections. Each population connects to 20 other populations, which are chosen randomly. There are thus 250 connections per neuron and projection. The way we scale up is by adding more populations. We still connect to 20 other populations.

#### To run MAM:

You need to have installed NEST *with* Python in order to run the MAM benchmark.

You also need to download *nested_dict* and *dicthash*:

```bash
pip install nested_dict --user
pip install dicthash --user
```

#### To run 4x4:

To install the 4x4 mesocircuit model, run

```bash
module load cray-python/2.7.15.1
module load PyExtensions/2.7.15.1-CrayGNU-18.08
module load h5py/2.8.0-CrayGNU-18.08-python2-parallel
pip install NeuroTools
```

Then, go to `nest-benchmarks/BenchModels/4x4mm2LFP/` and run

```bash
> python setup.py install --user
```

**Specifications:**
All of the benchmarks are run with the following number of nodes, VPs, scales, threads unless otherwise specified:

| Nodes | VPs  | Scale - hpc | Scale - hpc_syn | Scale - hpc_rule | Scale - pop | Scale - MAM | Scale - 4x4 | Threads |
| :---: | :--: | :---------: | :-------------: | :--------------: | :---------: | :---------: | :---------: | :-----: |
|   1   |  36  |     20      |      44.4       |        5         |      5      |   0.0625    |   0.0625    |    6    |
|   2   |  72  |     40      |      88.8       |        10        |     10      |    0.125    |    0.125    |    6    |
|   4   | 144  |     80      |      177.6      |        20        |     20      |    0.25     |    0.25     |    6    |
|   8   | 288  |     160     |      355.2      |        40        |     40      |     0.5     |     0.5     |    6    |
|  16   | 576  |     320     |      710.4      |        80        |     80      |     1.0     |     1.0     |    6    |
|  32   | 1152 |     640     |     1420.8      |       160        |     160     |     2.0     |     2.0     |    6    |



### Analyze Benchmarks

When the benchmark is finished running, tell JUBE to analyze the results

```bash
jube analyse <outpath> -i id
```

Then get the results

```bash
jube result <outpath> -i id
```

If you want to save the results in a csv-file, do

```bash
jube result <outpath> -i id > <path-to-directory-you-want>/result-name.csv
```

`result-name.csv` can be an existing file, or new one. If it already exists, this will overwrite the content in the file.

**Visualization**

[WIP]














