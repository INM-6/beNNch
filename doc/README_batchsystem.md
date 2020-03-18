
## Batch System config

    egrep -ho '#[A-Z_]+#' *.batch* | sort | uniq -c | sort -g

   Batch-system templates currently replace the following variables enclosed in '#':

### Allocation parameters

   `PARTITION`
   : batch partition to run on

   `WALLTIME`
   : maximum time of the allocation

   `MEMORY`
   : reserved maximum memory per node of the allocation

   `CONSTRAINT`
   : ???


### Job Geometry

   `TASKS`
   : number of processes started in the job

   `CPUS_PER_TASK`
   : number of compute-elements (cores) to reserve per process.

   `TASKS_PER_CORE`
   : number of processes per core

   `NODES`
   : number of nodes requested for the allocation

   `TASKS_PER_NODE`
   : number of processes per node


### Job Configuration

   `JOB_NAME`
   : the name of the job

   `JOBSCRIPT`
   : This is the actual script that runs the job. It should be agnostic of the calling batch system.

   `OUTPATH` and `ERRPATH`
   : path to stderr/stdout output file, may contain special characters that are interpreted by the batch system (e.g. Job id)

   `READY`
   : name of ready-file for JUBE to detect completion of the job.


### Batch system specific config

   Usage of these variables should be avoided as far as possible.

   `MAIL_ADDRESS`
   : notification mail address

   `MAIL_MODE`
   : notification mode


