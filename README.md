
# Requirements

  To run this workflow, you need to

  * have `conda` available,
  * have a `Snakemake` >=5.0 installed in your current environment, and
  * `$PYTHONPATH` must not be set.

  The workflow promises

  * to NOT modify the outside environment (e.g. by installing packages), and
  * to take care of all required dependencies internally.


# Installing

   ```bash
   # wget https://...miniconda3.sh
   # bash miniconda3.sh
   # conda install pip
   # which pip
   .../miniconda3/bin/pip
   # pip install snakemake
   # snakemake --version
   5.4.4
   ```

# Framework Description

   [...]

   When debugging, you can load commands to quickly jump into and out of a job
   result by sourcing `common/aliasses`

    ```bash
    source common/aliasses
    ```

   This will give you the commands `gostep` to jump into the workdir and see
   all `stderr` files, and `leavestep` to jump back to the base directory of
   the workflow.


   Look at <http://www.fz-juelich.de/ias/jsc/EN/Expertise/Support/Software/JUBE/JUBE2/jube-documentation_node.html>

   [Online Documentation](https://apps.fz-juelich.de/jsc/jube/jube2/docu/index.html)

