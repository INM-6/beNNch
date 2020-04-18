
# jube-environment-loading

   Environment setup depends on many external factors and is flexibly handled
   by the jube configuration. The environment can be set up by different means
   like `conda` or `module load`, etc. and the availability and package names
   are location dependant. To handle this in a flexible way the
   `simulator_params` in this example encode the various dependencies of which
   the benchmark selects the right one. Using the `init_with=` mechanism, a
   benchmark can easily specialize a parameter or define a scan (aka, template
   parameter).

    $ jube run envtest.xml
    $ cat envtest/000000/000000_run/work/stdout

   When changing the `machine_name` parameter, or the machine specific
   `envloader` and re-running:

    $ jube run envtest.xml
    $ find envtest/ -name stdout | sort | xargs head -n 100

   Where the second line is just an automated way of showing all `stdout` files
   produced so far. From time to time cleanup the results with

    $ rm -rf envtest/

