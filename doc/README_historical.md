
# NEST Benchmarking Framework

   This repository contains the config files used for running NEST-related
   benchmarks on various hosts. It uses the [Jülich Benchmarking
   Environment](https://apps.fz-juelich.de/jsc/jube/jube2/docu/index.html) to
   construct the parameter space to explore and build/run/analyse the different
   configurations on potentially big machines.

   The first set of benchmarks was developed at the JSC Performance Tunathon in
   February 2019.

# How to run

  First create a local clone of the relevant source code that contains all
  commits that you want to benchmark. This clone will be used by all benchmarks
  as source repostiroy for checkouts of specific revisions to benchmark.

    git clone git@github.com:INM-6/mini-nest

  This repository is in the same directory as the benchmark configuration in
  `test_bench.xml`. Do a benchmark run with:

    jube run test_bench.xml

  This will take some time. The data will be written into a newly created
  directory defined in the xml file (currently `test_bench/`). When the run is
  finished, the analysis can be performed and results displayed. Jube will tell
  you which commands to run (read the output!), e.g.:

    jube result test_bench/ -a --id 1

  For displaying results from several runs and maybe update the result table
  columns:

    jube result test_bench/ -u test_bench.xml -a --id 18 19


# Other repositories with Benchmark activities

  * ACA: https://jugit.fz-juelich.de/aca_requirements_validation_and_benchmarking/development-of-test-and-science-cases/tree/master/testcases
    * (ACA related) https://github.com/INM-6/spiking-htm
  * ICIE: https://wiki.humanbrainproject.eu/bin/view/Collabs/hbp-benchmark-suite-for-technology-trans/
  * Multi-Area Model: https://github.com/INM-6/mam_benches
  * Microcircuit: https://github.com/INM-6/microcircuit_jube_benches
  
  
# More information

   JUBE Documentation
   : <https://apps.fz-juelich.de/jsc/jube/jube2/docu/index.html>

