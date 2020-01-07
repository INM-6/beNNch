#!/usr/bin/env python
# encoding: utf8
from populations import CreatePopulations
from conn_rules import ConnectAll
import nest
import sys
import time
import logging
log = logging.getLogger()
logging.basicConfig(level=logging.DEBUG, filename='mylog.log', filemode='w')

args = sys.argv

if len(args) > 1:
    scale = int(args[1])
    totVPs = int(args[2])
else:
    scale = 0.1
    totVPs = 4

print("scale: ", scale)
print("totVPs: ", totVPs)

if __name__ == '__main__':
    nest.SetKernelStatus({'total_num_virtual_procs': totVPs})
    pops = CreatePopulations(scale)
    
    ConnectAll(pops)
    
    # Init time and memory
    tic = time.time()
    nest.Prepare()
    nest.Run(10.)

    InitializationTime = time.time() - tic

    print('{} # init_time'.format(InitializationTime))
    print('{} # virt_mem_after_init'.format(memory_thisjob()))

    print('{} # virt_mem_after_sim'.format(nest.ll_api.sli_func('memory_thisjob'))) # No simulation, just here for consistency
