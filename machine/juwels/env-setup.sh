#
# ENVIRONMENT SETUP
#
set -euo pipefail
module purge
# module load GCC CMake ParaStationMPI Python SciPy-Stack GSL git
module load #MODULES#
module list
set -x
export TIME='CMD: %C\nEXIT: %x\n\nMEM_UNSHARED: %DKiB\nMEM_STACK: %pKiB\nMEM_TOTAL: %KKiB\nMEM_RES: %MKiB\nMEM_SHARED: %XKiB\nMEM_RES_SET: %tKiB\nMEM_PAGEFAULTS: %F\nMEM_PAGERECOVERS: %R\n\nCPU_USER: %Us\nCPU_SYS: %Ss\nCPU_REAL: %es\nCPU_PCT: %P\n\nFS_INP: %I\nFS_OUT: %O\nSOCK_SEND: %smessages\nSOCK_RECV: %rmessages\nSWITCHES: %w\nSIGNALS: %k';
source compile/bin/nest_vars.sh
