#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Hybrid LFP scheme example script, applying the methodology with a model
implementation similar to:

Nicolas Brunel. "Dynamics of Sparsely Connected Networks of Excitatory and
Inhibitory Spiking Neurons". J Comput Neurosci, May 2000, Volume 8,
Issue 3, pp 183-208

But the network is implemented with spatial connectivity, i.e., the neurons
are assigned positions and distance-dependent connectivity in terms of
cell-cell connectivity and transmission delays.

Synopsis of the main simulation procedure:
1. Loading of parameterset
    a. network parameters
    b. parameters for hybrid scheme
2. Set up file destinations for different simulation output
3. network simulation
    a. execute network simulation using NEST (www.nest-initiative.org)
    b. merge network output (spikes, currents, voltages)
4. Create a object-representation that uses sqlite3 of all the spiking output 
5. Iterate over post-synaptic populations:
    a. Create Population object with appropriate parameters for
       each specific population
    b. Run all computations for populations
    c. Postprocess simulation output of all cells in population
6. Postprocess all cell- and population-specific output data
7. Create a tarball for all non-redundant simulation output

The full simulation can be evoked by issuing a mpirun call, such as
mpirun -np 4 python example_brunel.py

Not recommended, but running it serially should also work, e.g., calling
python example_brunel.py


Given the size of the network and demands for the multi-compartment LFP-
predictions using the present scheme, running the model on nothing but a large-
scale compute facility is strongly discouraged.


Copyright (C) 2014-2018, 4x4mm2LFP model authors.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

'''

from __future__ import division
import os
if 'DISPLAY' not in os.environ:
    import matplotlib
    matplotlib.use('Agg')
import sys
import numpy as np
import matplotlib.pyplot as plt
from  matplotlib.animation import FuncAnimation
from time import time
from hybridLFPy import PostProcess, setup_file_dest
from mesocircuit_analysis import CachedTopoNetwork, TopoPopulation, helpers
from mesocircuit_LFP_parameters import get_LFP_parameters
import parameterspace_control as psc
import h5py
import neuron
from mpi4py import MPI


#################################################
# matplotlib settings                           #
#################################################
plt.close('all')
plt.rcParams.update({'figure.figsize': [10.0, 8.0]})


#set some seed values
SEED = 12345678
SIMULATIONSEED = 12345678
np.random.seed(SEED)


#################################################
# Initialization of MPI stuff                   #
#################################################
COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()


#if True, execute full model. If False, do only the plotting.
#Simulation results must exist.
PROPERRUN = True


#check if mod file for synapse model specified in alphaisyn.mod is loaded
if not hasattr(neuron.h, 'ExpSynI'):
    if RANK == 0:
        os.system('nrnivmodl')
    COMM.Barrier()
    neuron.load_mechanisms('.')


#################################################################################
### HYBRID SCHEME PARAMETERS
#################################################################################

#Fetch parameter set file and corresponding ID
#if no additional argvs are supplied, fall back to base parameters
parameter_set_file = sys.argv[-1]
if parameter_set_file.endswith('python') or parameter_set_file.endswith('.py'):
    ps_id = helpers.get_unique_id(psc.base_parameters)
    parameter_set_file = os.path.join('parameters', ps_id + '.sli')
else:
    ps_id = os.path.split(sys.argv[-1])[-1][:-4]

#get corresponding parameterSet dictionary
PS, PSET = get_LFP_parameters(parameter_set_file, ps_id)


################################################################################
# MAIN simulation procedure                                                    #
################################################################################

if PROPERRUN:
    #set up the file destination, removing old results by default
    setup_file_dest(PS, clearDestination=False)

#wait for operation to finish
COMM.Barrier()

if RANK == 0:
    simstats = file(os.path.join(PS.savefolder, 'simstats.dat'), 'w')


#tic toc
tic = time()
ticc = tic

#Create an object representation containing the spiking activity of the network
#simulation output that uses sqlite3. 
networkSim = CachedTopoNetwork(**PS.network_params)

tocc = time()
if RANK == 0:
    simstats.write('CachedNetwork {}\n'.format(tocc - ticc))
    
toc = time() - tic
print(('NEST simulation and gdf file processing done in  %.3f seconds' % toc))



####### Set up populations #####################################################

if PROPERRUN:
    #iterate over each cell type, and create populationulation object
    for i, y in enumerate(PS.y):
        # if y in ['b4','nb4', 'p5(L23)','p5(L56)', 'b5','nb5', 'p6(L4)','p6(L56)', 'b6','nb6']:
        # if y in ['p6(L56)', 'b6','nb6']:
        #create population:
        ticc = time()
        pop = TopoPopulation(
                cellParams = PS.cellParams[y],
                rand_rot_axis = PS.rand_rot_axis[y],
                simulationParams = PS.simulationParams,
                populationParams = PS.populationParams[y],
                y = y,
                layerBoundaries = PS.layerBoundaries,
                electrodeParams = PS.electrodeParams,
                savelist = PS.savelist,
                savefolder = PS.savefolder,
                calculateCSD = PS.calculateCSD,
                dt_output = PS.dt_output, 
                POPULATIONSEED = SIMULATIONSEED + i,
                X = PS.X,
                networkSim = networkSim,
                k_yXL = PS.k_yXL[y],
                synParams = PS.synParams[y],
                synDelayLoc = PS.synDelayLoc[y],
                synDelayScale = PS.synDelayScale[y],
                J_yX = PS.J_yX[y],
                tau_yX = PS.tau_yX[y],
                #TopoPopulation kwargs
                topology_connections = PS.topology_connections,
                #output file name(s)
                output_file = PS.output_file,
            )
        tocc = time()
        if RANK == 0:
            simstats.write('Population_{} {}\n'.format(y, tocc-ticc))
    
        #run population simulation and collect the data
        ticc = time()
        pop.run()
        tocc = time()
        if RANK == 0:
            simstats.write('run_{} {}\n'.format(y, tocc-ticc))

        ticc = time()
        pop.collect_data()
        tocc = time()
        if RANK == 0:
            simstats.write('collect_{} {}\n'.format(y, tocc-ticc))
    
        #object no longer needed
        del pop


####### Postprocess the simulation output ######################################

#reset seed, but output should be deterministic from now on
np.random.seed(SIMULATIONSEED)

if PROPERRUN:
    ticc = time()
    #do some postprocessing on the collected data, i.e., superposition
    #of population LFPs, CSDs etc
    postproc = PostProcess(y = PS.y,
                           dt_output = PS.dt_output,
                           savefolder = PS.savefolder,
                           mapping_Yy = PS.mapping_Yy,
                           savelist = PS.savelist,
                           cells_subfolder = os.path.split(PS.cells_path)[-1],
                           populations_subfolder = os.path.split(PS.populations_path)[-1],
                           figures_subfolder = os.path.split(PS.cells_path)[-1],
                           output_file = PS.output_file,
                           compound_file = PS.compound_file,
                           )
    
    #run through the procedure
    postproc.run()
    
    #create tar-archive with output for plotting, ssh-ing etc.
    postproc.create_tar_archive()
    
    tocc = time()
    if RANK == 0:
        simstats.write('postprocess {}\n'.format(tocc-ticc))
        simstats.close()
    
COMM.Barrier()

#tic toc
print(('Execution time: %.3f seconds' %  (time() - tic)))
