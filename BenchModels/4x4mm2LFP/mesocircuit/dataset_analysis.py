#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Analysis routine for simulated data.

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
import h5py
import os
from time import time, sleep
from mesocircuit_analysis import helpers, NetworkAnalysis
from parameterspace import base_parameters
import parameterspace_control as psc
from mpi4py import MPI
import copy
import scipy.sparse as sp


###################################
# Initialization of MPI stuff     #
###################################
COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()


def get_network_analysis_object(parameter_set_file, ps_id,
                                TRANSIENT=0.,
                                BINSIZE_TIME=0.5,
                                BINSIZE_AREA=0.1):
    '''
    Return the meso_analysis.NetworkAnalysis object corresponding to a
    parameterset hash id
    
    Implemented as a mean to avoid defining parameters such as TRANSIENT,
    BINSIZE_TIME, BINSIZE_AREA in more than one place.
    
    Parameters
    ----------
    parameter_set_file : str
        path to parameterset sli file
    ps_id : str
        unique id of parameterset
    
    Returns
    -------
    object : meso_analysis.analysis.NetworkAnalysis
    '''
    #set up the analysis class instance, it is now a daughter class of
    #NetworkAnalysisParams, and will pass kwargs to parent class
    analysis = NetworkAnalysis(
                 parameter_set_file=parameter_set_file,
                 output_raw_prefix=os.path.join(psc.output_raw_prefix, ps_id),
                 output_proc_prefix=os.path.join(psc.output_proc_prefix, ps_id),
                 TRANSIENT=TRANSIENT,
                 BINSIZE_TIME=BINSIZE_TIME,
                 BINSIZE_AREA=BINSIZE_AREA,
    )
    
    return analysis


def write_sparse_datasets_to_h5_X(analysis, X, filenames, datasets):
    '''
    Writes sparse datasets for population X to h5 file.

    Parameters
    ----------
    analysis : meso_analysis.analysis.NetworkAnalysis
    X : population name
    filenames : list of filenames
    datasets : list of datasets
    '''
    for fname, data in zip(filenames, datasets):
        while True:
            try:
                fpath = os.path.join(analysis.output_proc_prefix, fname + '_{0}.h5'.format(X))
                f = h5py.File(fpath, 'w')
                helpers.dump_sparse_to_h5(X, f=f, data=data,
                                          compression='gzip', compression_opts=2)
                f.flush()
                f.close()
                break

            except: # IOError as e:
                print("Expected bad hdf5 behaviour:", sys.exc_info(), fpath)
                sleep(1.)
    return


def write_datasets_to_h5_X(analysis, X, filenames, datasets, dtypes):
    '''
    Writes non-sparse datasets for population X to h5 file.

    Parameters
    ----------
    analysis : meso_analysis.analysis.NetworkAnalysis
    X : population name
    filenames : list of filenames
    datasets : list of datasets
    dtypes : list of datatypes of data in datasets
    '''
    for fname, data, dtype in zip(filenames, datasets, dtypes):
        fpath = os.path.join(analysis.output_proc_prefix, fname + '_{0}.h5'.format(X))
        while True and data is not None:
            try:
                f = h5py.File(fpath, 'w')

                if np.isscalar(data):
                    f.create_dataset(X, data=data)
                else:
                    if data.shape != (0,): # needed for python 2.6
                        f.create_dataset(X, data=data,
                                         dtype=dtype, compression='gzip',
                                         compression_opts=2,
                                         chunks=True,
                                         shape=data.shape)
                f.flush()
                f.close()
                break

            except IOError as e:
                print("Expected bad hdf5 behaviour:", sys.exc_info(), fpath)
                sleep(np.random.rand())
    return


def combine_files_for_populations(analysis, filenames):
    '''
    Combines h5 files for different populations.

    Parameters
    ----------
    analysis : meso_analysis.analysis.NetworkAnalysis
    filenames : list of filenames without ending _X
    '''
    for fname in filenames:
        fpath1 = os.path.join(analysis.output_proc_prefix, fname + '.h5')
        try:
            f1 = h5py.File(fpath1, 'w')
            for X in analysis.X:
                fpath0 = os.path.join(analysis.output_proc_prefix, fname + '_{0}.h5'.format(X))
                f0 = h5py.File(fpath0, 'r')
                if X in f0.keys(): # quick check if key exists
                    f1.copy(f0[X], X)
                f0.close()
                os.system('rm {0}'.format(fpath0))
            f1.close()
        except IOError as ioe:
            print 'error', fpath0
            f1.close()
            os.remove(fpath1)
    return


if __name__ == '__main__':
    #tic toc
    tic = time()

    
    #if no additional argvs are supplied, fall back to base parameters
    parameter_set_file = sys.argv[-1]
    if parameter_set_file.endswith('python') or parameter_set_file.endswith('.py'):
        ps_id = helpers.get_unique_id(base_parameters)
        parameter_set_file = os.path.join('parameters', ps_id + '.sli')
    else:
        ps_id = os.path.split(sys.argv[-1])[-1][:-4]
    
    analysis = get_network_analysis_object(parameter_set_file, ps_id)
    
    #run through the main nest output files, collect GDFs, spikes and positions
    analysis.run()

    # firing rates
    analysis.print_firing_rates()
    
    #compute grid patch specific histograms of neuron locations, spike counts,
    #neuron spike rates, time-downsampled rates


    #iterate over all network populations
    print('Computing results of pop: ')
    #argsort as simple attempt at loadbalancing, also across compute nodes:
    for i, j in enumerate(np.argsort(analysis.N_X)[::-1]):
        X = analysis.X[j]
        if i % SIZE == RANK:
            print('population {0} on RANK {1}'.format(X, RANK))
            
            ####################################################################
            # PREPROCESS DATA
            ####################################################################
            #get the spike output of neurons with corrected GIDs
            spikes = helpers.read_gdf(os.path.join(analysis.output_proc_prefix,
                                                analysis.spike_detector_label +
                                                X + '.gdf'))
            
            #positions
            binned_positions = analysis.compute_position_hist(analysis.positions_corrected[X])
            
    
            #position binned spike trains
            sptrains = analysis.compute_time_binned_sptrains(X, spikes, analysis.time_bins, dtype=np.uint8)
            binned_sptrains = analysis.compute_pos_binned_sptrains(analysis.positions_corrected[X],
                                                                   sptrains,
                                                                   dtype=np.uint16)
            
            #spike counts and rates per spatial bin
            binned_spcounts = np.asarray(binned_sptrains.sum(axis=1).astype(int).reshape((analysis.pos_bins.size-1, -1)))
            binned_mean_sprates = analysis.compute_binned_rates(binned_spcounts, binned_positions)

    

            #create some time resampled spike and rate histograms
            sptrains_rs = analysis.compute_time_binned_sptrains(X, spikes,
                                                                analysis.time_bins_rs,
                                                                dtype=np.uint8)
            binned_sptrains_rs = analysis.compute_pos_binned_sptrains(analysis.positions_corrected[X],
                                                                      sptrains_rs,
                                                                      dtype=np.uint16)
            #spatially binned resampled spike rates
            binned_sprates_rs = analysis.compute_binned_sprates(binned_sptrains_rs,
                                                                binned_positions,
                                                                analysis.time_bins_rs)
            
            
            #position sorting arrays
            pos_sorting_arrays = analysis.get_pos_sorting_array(analysis.positions_corrected[X])


            ####################################################################
            # WRITE STUFF
            ####################################################################
            
            #condense syntax for writing
            filenames = [
                'all_sptrains',
                'all_sptrains_rs',
                'all_binned_sptrains',
                'all_binned_sptrains_rs',
                'all_binned_sprates_rs',
                ]
            datasets = [
                sptrains,
                sptrains_rs,
                binned_sptrains,
                binned_sptrains_rs,
                binned_sprates_rs,
                ]

            write_sparse_datasets_to_h5_X(analysis, X, filenames, datasets)


            filenames = [
                'all_binned_positions',
                'all_binned_spcounts',
                'all_pos_sorting_arrays',
                'all_binned_mean_sprates',
            ]

            datasets = [
                binned_positions,
                binned_spcounts,
                pos_sorting_arrays,
                binned_mean_sprates,
            ]

            dtypes = [int for _ in range(3)] + [float]

            write_datasets_to_h5_X(analysis, X, filenames, datasets, dtypes)

    COMM.Barrier()


    #put certain output data as generated on different parallel threads in
    #single files for further analysis.
    if RANK == 0:
        filenames = [
            'all_sptrains',
            'all_sptrains_rs',
            'all_binned_positions',
            'all_binned_spcounts',
            'all_binned_sptrains',
            'all_binned_sptrains_rs',
            'all_pos_sorting_arrays',
            'all_binned_mean_sprates',
            'all_binned_sprates_rs',
        ]

        combine_files_for_populations(analysis, filenames)

    COMM.Barrier()

    if RANK == 0:
        toc = time()-tic
        print('analysed in {0} seconds'.format(toc))

        #write analysis_info.dat file with analysis-, plotting- and animation times
        f = file(os.path.join(analysis.output_proc_prefix, 'analysis_info.dat'), 'w')
        f.write('analysis: %.2f s\n' % toc)
        f.close()

    COMM.Barrier()
