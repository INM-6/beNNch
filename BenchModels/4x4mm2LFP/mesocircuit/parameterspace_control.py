#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Creation of parameter files and job scripts for mesocircuit model.

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
import sys
import pickle as pickle
from NeuroTools import parameters as ps
from mesocircuit_analysis import helpers
from parameterspace import base_parameters, ParameterSpaces
import data_areas_species
import derive_params as derive
from mean_conn_prob_area_size import get_mean_delay_circle

args = sys.argv
print(args)
scale = float(args[1])
totVPs = int(args[2])
simtime = float(args[3])
transient = float(args[4])

###########################################################
# set up jobscripts for different machines                #
###########################################################

if False: #'HOSTNAME' in os.environ:
#-------------------------------------------------------------------------------
    if os.environ['HOSTNAME'].rfind('blaustein') >= 0 or os.environ['HOSTNAME'
                                                                      ].rfind('jr') >= 0:
        jobscript_skeleton_network = '''#!/bin/bash
#SBATCH --job-name %s
#SBATCH --time %s
#SBATCH -e %s
#SBATCH -o %s
#SBATCH -N %i
#SBATCH --cpus-per-task %i
#SBATCH --exclusive
cd %s
%s nest %s mesocircuit.sli
'''
        jobscript_skeleton_analysis = '''#!/bin/bash
#SBATCH --job-name %s
#SBATCH --time %s
#SBATCH -e %s
#SBATCH -o %s
#SBATCH -N %i
#SBATCH --cpus-per-task %i
#SBATCH --exclusive
unset DISPLAY
cd %s
%s python dataset_analysis.py %s
'''
        jobscript_skeleton_LFP = '''#!/bin/bash
#SBATCH --job-name %s
#SBATCH --time %s
#SBATCH -e %s
#SBATCH -o %s
#SBATCH -N %i
#SBATCH --ntasks %i
#SBATCH --exclusive
unset DISPLAY
cd %s
%s python %s %s
'''

#-------------------------------------------------------------------------------
    elif os.environ['HOSTNAME'].rfind('local') >= 0:
        jobscript_skeleton_network = '''#!/bin/sh
#PBS -o %s
#PBS -e %s
#PBS -lnodes=%i:ppn=%i
#PBS -lwalltime=%s
cd $PBS_O_WORKDIR
mpirun -np %i -bynode -bind-to-core -cpus-per-proc %i nest %s mesocircuit.sli
wait
'''
        jobscript_skeleton_analysis = '''#!/bin/sh
#PBS -o %s
#PBS -e %s
#PBS -lnodes=%i:ppn=%i
#PBS -lwalltime=%s
cd $PBS_O_WORKDIR
mpirun -pernode python dataset_analysis.py %s
'''
        jobscript_skeleton_LFP = '''#!/bin/sh
#PBS -o %s
#PBS -e %s
#PBS -lnodes=%i:ppn=%i
#PBS -lwalltime=%s
cd $PBS_O_WORKDIR
mpirun -np %i python %s %s
wait
'''
    else:
        raise Exception('HOSTNAME {0} not recognized, won\'t write jobscripts'.format(os.environ['HOSTNAME']))

#-------------------------------------------------------------------------------
else: #write nothing!
    jobscript_skeleton_network = None
    jobscript_skeleton_analysis = None
    jobscript_skeleton_LFP = None
#-------------------------------------------------------------------------------
    



###################################################
### Global  simulation parameters               ###
###################################################

#check where I am
if 'HOSTNAME' in os.environ:
    if os.environ['HOSTNAME'].rfind('blaustein') >= 0 or \
            os.environ['HOSTNAME'].rfind('jr') >= 0 or \
            os.environ['HOSTNAME'].rfind('local') >= 0:
        I_AM_SUPERCOMPUTER = True
    else:
        I_AM_SUPERCOMPUTER = False
else:
    I_AM_SUPERCOMPUTER = False

#flag for whether or not jobs should be rerun (on basis of existance of some files)
#NOTE: so far only tested to work on Stallo
FORCE_RERUN = True

# slurm jobscript parameters
if I_AM_SUPERCOMPUTER:
    if os.environ['HOSTNAME'].rfind('local') >= 0:
        lnodes = 16  #number of compute nodes for network sim
        lnodes_analysis = 4 #number of nodes for analysis
        lnodes_LFP = 15 #number of nodes for LFP calculations
        ppn = 16    #number of cores per node
        cpus_per_proc = ppn #CPUs per MPI process, use this to optimize memory
        mem_per_proc = 2000  #mb, memory per CPU core
        num_rec_procs = 1 #number of recording mpi processes, use GSD if > 0
        output_raw_prefix = os.path.join('/', 'global', 'work',
                                         os.environ['USER'], 'mesocircuit',
                                         'output_raw')
        output_proc_prefix = os.path.join('/', 'global', 'work',
                                          os.environ['USER'], 'mesocircuit',
                                          'output_processed')
        jobsubmit = 'qsub'
        n_vp = lnodes*ppn*8
        lwalltime_network = '0:30:00'
        lwalltime_LFP = '6:00:00'
    elif os.environ['HOSTNAME'].rfind('blaustein') >= 0:
        lnodes = 10         #number of node boards
        lnodes_analysis = 2 #number of nodes for analysis
        lnodes_LFP = 16 #number of nodes for LFP calculations
        ppn = 48            #(virtual) processors per node
        mem_per_proc = 64000 #MB memory per process, not used
        num_rec_procs = 1   #number of recording processes
        output_raw_prefix = 'output_raw'
        output_proc_prefix = 'output_processed'
        jobsubmit = 'sbatch'
        n_vp = lnodes*ppn
        lwalltime_network = '0:30:00'
        lwalltime_analysis = '0:30:00'
        lwalltime_LFP = '6:00:00'
    elif os.environ['HOSTNAME'].rfind('jr') >= 0:
        lnodes = 48         #number of node boards
        lnodes_analysis = 2 #number of nodes for analysis
        lnodes_LFP = 96 #number of nodes for LFP calculations
        ppn = 48            # number of (virtual) processes per node
        mem_per_proc = 16000 #MB memory per process, not used
        num_rec_procs = 1   #number of recording processes
        output_raw_prefix = os.path.join('/', 'work', 'jinb33',
                                         os.environ['USER'], 'mesocircuit',
                                         'output_raw')
        output_proc_prefix = os.path.join('/', 'work', 'jinb33',
                                         os.environ['USER'], 'mesocircuit',
                                         'output_processed')
        jobsubmit = 'sbatch'
        n_vp = lnodes*ppn # number of virtual prcocesses, assume multithreading
        lwalltime_network = [20, 200, 200, 20] #[create, connect, sim per s, write], measured on jureca
        lwalltime_analysis = [400, 80, 70, 500] #[analysis per s, utah per s, plotting per s, animation], measured on jureca
        lwalltime_LFP = [25, 1000, 15000, 1000, 600] #[cache per s, pop setup, sim per s, collect per s, postprocess], jureca
    else:
        #presumably running on local computer
        lnodes = 1  #number of compute nodes
        lnodes_analysis = 1 #number of nodes for analysis
        lnodes_LFP = 1 #number of nodes for LFP calculations
        ppn = 4    #number of cores per node
        cpus_per_proc = 1 #CPUs per MPI process, use this to optimize memory
        mem_per_proc = 2000  #mb, memory per CPU core
        num_rec_procs = 1 #number of recording mpi processes, use GSD if > 0
        output_raw_prefix = 'output_raw'
        output_proc_prefix = 'output_processed'
        lwalltime_network = '0:30:00'
        lwalltime_analysis = '0:30:00'
        lwalltime_LFP = '6:00:00'

else:
    lnodes = 1  #number of compute nodes
    lnodes_analysis = 4 #number of nodes for analysis
    lnodes_LFP = 1 #number of nodes for LFP calculations
    ppn = 4    #number of cores per node
    cpus_per_proc = 1 #CPUs per MPI process, use this to optimize memory
    mem_per_proc = 1000  #mb, memory per CPU core
    num_rec_procs = 0 #number of recording mpi processes, use GSD if > 0
    output_raw_prefix = 'output_raw'
    output_proc_prefix = 'output_processed'
    n_vp = 4
    lwalltime_network = '0:30:00'
    lwalltime_analysis = '0:30:00'
    lwalltime_LFP = '6:00:00'


#estimated running time including python postprocessing and plotting
# lwalltime_network = [200, 600, 150, 100] #[create, connect, sim per s, write], measured on Juqueen     
# lwalltime_network = '0:30:00'    
# lwalltime_analysis = '0:30:00'    
# lwalltime_LFP = '6:00:00'    
lwalltime_LFP_proxy = '2:30:00' #3900 s measured on jureca, fixed time  (7200 w. periodic boundaries)  


#number of virtual processes, must be even dividable by
#MPI SIZE, now also assigning two threads per core for small performance boost
#n_vp = 2048 #this is what we use on Juqueen, 128*16
#n_vp = lnodes*ppn*2


# File destinations
parameterset_dest = os.path.dirname(__file__) +  '/parameters/' + str(totVPs) + '-' + str(scale) + '/'
jobscript_dest = 'jobscripts'
log_dir = 'log'


def write_jobscripts(jobscript_skeleton_network,
                     jobscript_skeleton_analysis,
                     jobscript_skeleton_LFP):
    '''Deal with the parameterSpace object, doing file i/o, create jobs etc.'''
    #set up some file destination folders
    if not os.path.isdir(parameterset_dest):
        os.mkdir(parameterset_dest)

    if not os.path.isdir(jobscript_dest):
        os.mkdir(jobscript_dest)

        
    #dump base_parameters as txt file
    base_parameters.save('base_parameters.txt')
    
    #container
    PSdict = {}

    #iterate over different parameterspaces in dict
    for key, PS in list(ParameterSpaces.items()):

        #save ParameterSpace object using built-in method
        PS.save('PS_{0}.txt'.format(key))
    
    
        for paramset in PS.iter_inner():
            #unique id for each parameter set, constructed from the parameset dict
            #converted to a sorted list of tuples
            ps_id = helpers.get_unique_id(paramset)
            
            if ps_id in list(PSdict.keys()):
                print('skipping {0}, already in job list'.format(ps_id))
                pass
            else:
                print('Hash: {0}'.format(ps_id))
                #generated sli code goes here
                parameter_set_file = os.path.join(parameterset_dest, 
                                                  '%s.sli' % ps_id)
        
                #unique output path for each parameter set
                output_path = os.path.join(output_raw_prefix, ps_id)
                
                if not os.path.isdir(output_raw_prefix):
                    os.mkdir(output_raw_prefix)
                if not os.path.isdir(output_path):
                    os.mkdir(output_path)
                # just the prefix folder
                if not os.path.isdir(output_proc_prefix):
                    os.mkdir(output_proc_prefix)
                
            
                #put output_path into dictionary, as we now have a unique ID of
                #though this will not affect the parameter space object PS
                paramset = paramset.copy()
                paramset.update({
                    'n_vp' : n_vp,
                    'ps_id' : ps_id,
                    'num_rec_procs' : num_rec_procs,
                    'output_path' : output_path,
                    })

                paramset.update({
                    'density': scale,
                    't_sim': simtime,
                    'transient': transient,
                    'n_vp': totVPs
                    })

                paramset['t_sim'] = simtime

                ################################################################
                # Block for parameters derived from other parameters that may
                # be a parameterRange in a ParameterSpace, e.g., g
                ################################################################
                paramset = paramset.copy()

                # sets parameters: full_num_synapses, full_num_synapses_th,
                # weight_mod_facts, K_bg
                paramset.update(derived_params = derive.get_all_derived_params(paramset))
                # strip paramset of no-longer used parameters
                #paramset.pop('conn_prob_modifications')
                #paramset.pop('weight_modifications')


                # model-neuron parameter dictionary
                if paramset['neuron_model'] == '/iaf_psc_exp':
                    # add neuron parameters to dictionary
                    # (currently, tau_syn_ex and tau_syn_in are defined
                    # additionally in parameterspace and have to be copied in)
                    paramset['model_params'] = {}
                    params = ['tau_m', 'tau_syn_ex', 'tau_syn_in', 't_ref',
                              'E_L', 'V_th', 'C_m', 'V_reset']
                    for p in params:
                        paramset['model_params'].update({'/'+p : paramset[p]})
                        # mesocircuit_meanfield_predictions.py currently needs single neuron parameters
                        #paramset.pop(p)
                else:
                    raise NotImplementedError, 'model neuron {} not implemented'.format(paramset['neuron_model'])


                # get parameters for a specific spatial profile
                if paramset['spatial_profile'] != 'default':
                    for key in data_areas_species.spatial_profiles.keys():
                        paramset.update({key: \
                        data_areas_species.spatial_profiles[key][paramset['spatial_profile']]})

                # overwrite normal delay with mean of linear delay on disk of 1mm^2
                # and a decay parameter of 0.3 mm
                if paramset['use_mean_delays'] == True:
                    delays = []
                    for i,d0 in enumerate(paramset['min_delays']):
                        dmean = get_mean_delay_circle(area=1.,
                                                      d0=d0,
                                                      v=paramset['prop_speeds'][i],
                                                      sigma=0.3)
                        delays.append(dmean)
                    paramset['delays'] = delays


                #### END of BLOCK ######
                
                #put parameter set into dict using ps_id as key
                PSdict[ps_id] = paramset

                # update with parameters in 'derived_params' and then remove this
                # item
                PSdict[ps_id].update(PSdict[ps_id]['derived_params'])
                PSdict[ps_id].pop('derived_params')
            
                #put parameter set into executable sli file
                helpers.dump_sli_params(parameter_set_file, paramset)

                # write parameters to text file
                helpers.write_params(os.path.join(parameterset_dest, 
                                                  '%s.txt' % ps_id),
                                     PSdict[ps_id])
                
                #write using ParemeterSet native format, as we may need to load
                #exactly what was defined in the first place.
                ps.ParameterSet(PSdict[ps_id]).save(url=os.path.join(parameterset_dest,
                                                                     '{0}.pset'.format(ps_id)))
                
                #specify where to save output and errors
                if not os.path.isdir(log_dir):
                    os.mkdir(log_dir)
                output_network = os.path.join(os.getcwd(), log_dir, ps_id+'_network.txt')
                output_analysis = os.path.join(os.getcwd(), log_dir, ps_id+'_analysis.txt')
                output_LFP = os.path.join(os.getcwd(), log_dir, ps_id+'_LFP.txt')

    
    #dump dictionary of unique parametersets
    f = open('PSdict.pickle', 'w')
    pickle.dump(PSdict, f)
    f.close()


    print('created {0} unique jobs'.format(len(list(PSdict.keys()))))


if __name__ == '__main__':
    write_jobscripts(jobscript_skeleton_network,
                     jobscript_skeleton_analysis,
                     jobscript_skeleton_LFP)
