#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Main analysis functions.

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

#import sys
import numpy as np
import scipy.sparse as sp
import glob
import h5py
import json
import os
import copy
import math
from . import helpers
from mpi4py import MPI


###################################
# Initialization of MPI stuff     #
###################################
COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()


##################################
# Class definitions              #
##################################

class NetworkAnalysisParams(object):
    '''dummy class gathering some necessary simulation parameters'''
    def __init__(self,
                 parameter_set_file,
                 output_raw_prefix='output_raw',
                 output_proc_prefix='output_processed',
                 TRANSIENT=0.,
                 BINSIZE_TIME=1.,
                 BINSIZE_AREA=0.1,
                 clean_proc_prefix=False,
                ):
        '''
        initialization of class NetworkAnalysisParams

        Most needed parameter values are extracted from sli files

        Arguments
        ---------
        parameter_set_file : str
            /path/to/parameters.sli file
        output_raw_prefix :
            /path/to/raw_data
        output_proc_prefix : str
            /path/to/processed_data
        TRANSIENT : float
            starting point of analysis (to avoid start trans)
        BINSIZE_TIME : float
            binsize when temporally binning e.g., spikes
        BINSIZE_AREA : float
            binsize when spatially binning e.g., positions
        clean_proc_prefix : bool
            delete all existing files in /path/to/processed_data
        '''
        #set class attributes
        self.parameter_set_file = parameter_set_file  #read in from command line
        self.output_raw_prefix = output_raw_prefix
        self.output_proc_prefix = output_proc_prefix
        self.TRANSIENT = TRANSIENT
        self.BINSIZE_TIME = BINSIZE_TIME
        self.BINSIZE_AREA = BINSIZE_AREA
        self.clean_proc_prefix = clean_proc_prefix


        # Various parameters loaded from sli files
        self.record_fraction_neurons_spikes = helpers.get_sli_variable(
                    self.parameter_set_file, 'record_fraction_neurons_spikes')
        self.frac_rec_spikes = float(helpers.get_sli_variable(self.parameter_set_file,
                                                          'frac_rec_spikes'))
        self.n_rec_spikes = int(helpers.get_sli_variable(self.parameter_set_file,
                                                     'n_rec_spikes'))


        #time resolution of simulation
        self.dt = float(helpers.get_sli_variable(self.parameter_set_file, 'dt'))

        # simulation time
        self.t_sim = float(helpers.get_sli_variable(self.parameter_set_file,
                                                    't_sim'))

        # transient from NEST simulation, not to confuse with self.TRANSIENT,
        # the transient used just for the analysis
        self.transient = float(helpers.get_sli_variable(self.parameter_set_file,
                                                        'transient'))

        # get extent_length of the topology layer
        # positions are in [-0.5*extent_length, 0.5*extent_length]
        self.extent_length = float(helpers.get_sli_variable(self.parameter_set_file,
                                                    'extent_length'))

        # whether periodic boundary conditions were used
        self.pbc = helpers.get_sli_variable(self.parameter_set_file, 'pbc')


        #file name of raw population 
        self.GID_filename = helpers.get_sli_variable(self.parameter_set_file,
                                             '/GID_filename')


        #file name of GID positions
        self.position_filename = helpers.get_sli_variable(self.parameter_set_file,
                                                  '/position_filename_trunk')

        #file name of recorded connections
        self.connection_filename = helpers.get_sli_variable(self.parameter_set_file,
                                                    '/connection_filename_trunk')#[1:-1]


        #file prefix for spike detector
        self.spike_detector_label = helpers.get_sli_variable(self.parameter_set_file,
                                                '/spike_detector_label')

        #whether voltages have been saved
        self.save_voltages = helpers.get_sli_variable(self.parameter_set_file,
                                                'save_voltages')

        if self.save_voltages:
            #file prefix for voltmeter
            self.voltmeter_label = helpers.get_sli_variable(self.parameter_set_file,
                                                    '/voltmeter_label')

            #time interval for probing voltages
            self.dt_voltage = float(helpers.get_sli_variable(self.parameter_set_file,
                                                    'dt_voltage'))

        #presynaptic population names
        self.X = ['L23E', 'L23I', 'L4E', 'L4I',
                  'L5E', 'L5I', 'L6E', 'L6I', 'TC']


        #postsynaptic population names
        self.Y = self.X[:-1]


        #sorting axis, 'x', 'y' or None
        self.sorting_axis = 'x'


        #time bins for spike trains
        self.time_bins = np.arange(self.TRANSIENT / self.dt,
                            self.t_sim / self.dt) * self.dt
        self.time_bins_rs = np.arange(self.TRANSIENT / self.BINSIZE_TIME,
                            self.t_sim / self.BINSIZE_TIME) * self.BINSIZE_TIME


        #bins for the positions
        self.pos_bins = np.linspace(-self.extent_length / 2,
                                    self.extent_length / 2,
                                    int(self.extent_length /\
                                        self.BINSIZE_AREA + 1))
        
        

class NetworkAnalysis(NetworkAnalysisParams):
    '''main class NetworkAnalysis doing the analysis steps'''
    def __init__(self,
                 **kwargs
                 ):
        '''
        initialization of class NetworkAnalysis

        consider just doing this by inheritance of class NetworkAnalysisParams

        Arguments
        ---------
        kwargs : 
            see class NetworkAnalysisParams
        '''
        #initalize parent class
        NetworkAnalysisParams.__init__(self, **kwargs)


        # load GIDs
        GID_filename = open(os.path.join(self.output_raw_prefix,
                                         self.GID_filename), 'r')
        self.GIDs = []
        for l in GID_filename:
            a = l.split()
            self.GIDs.append([int(a[0]), int(a[1])])
        GID_filename.close()


        # population sizes
        self.N_X = [self.GIDs[i][1] - self.GIDs[i][0] + 1 # +1 added
                          for i in range(len(self.GIDs))]


        # discard population name 'TC' if no thalamic neurons are simulated
        self.X = self.X[:len(self.N_X)]


        # numbers of neurons for which spikes were recorded
        if self.record_fraction_neurons_spikes:
            self.rec_sizes = [int(self.N_X[i] * self.frac_rec_spikes)
                              for i in range(len(self.X))]
        else:
            self.rec_sizes = [self.n_rec_spikes] * len(self.X)

        # currently, spikes are recorded from all thalamic neurons
        self.rec_sizes[-1] = self.N_X[-1]

        if RANK == 0:
            print('\nGlobal uncorrected ids ([first id, last id]):')
            for GID in self.GIDs:
                print(GID)
            print('\n')
            
            print('Population sizes:')
            for X, s in zip(self.X, self.N_X):
                print('{0}:\t{1}'.format(X, s))
            print('\n')
            
            # total number of neurons in the simulation
            print('Total number of neurons:')
            print(sum(self.N_X))
            print('\n')
            

    
    def run(self):
        '''
        Default procedure for converting file output from nest
        '''
        # clear destination of processed nest output
        if RANK == 0:
            if os.path.isdir(self.output_proc_prefix):
                if self.clean_proc_prefix:
                    os.system('rm -rf %s' % self.output_proc_prefix)
                    os.mkdir(self.output_proc_prefix)
            else:
                os.mkdir(self.output_proc_prefix)
        #sync
        COMM.Barrier()
        
        # load spikes from gdf files, correct GIDs,
        # merge them in separate files, and store spike trains
        self.GIDs_corrected = self.get_GIDs()

        #combine spike files generated on individual threads into one file
        #per population
        self.merge_files(self.spike_detector_label)

        #get the positions, correct the GIDs and write to file:
        self.positions_corrected = self.get_positions()

        #get the recorded connections, correct the GIDs and write to file
        self.connections_corrected = self.get_connections()
        
        #get the spike data, fill in dict
        self.spike_detector_output = {}
        for i, X in enumerate(self.X):
            self.spike_detector_output[X] = helpers.read_gdf(os.path.join(
                                        self.output_proc_prefix,
                                        self.spike_detector_label + X + '.gdf'))

        if self.save_voltages:
            #combine voltage files generated on individual threads into one file
            #per population
            self.merge_files(self.voltmeter_label)

            #get the voltage data, fill in dict
            self.voltmeter_output = {}
            for i, X in enumerate(self.X):
                self.voltmeter_output[X] = helpers.read_gdf(os.path.join(
                    self.output_proc_prefix,
                    self.voltmeter_label + X + '.gdf'))
        


    def get_GIDs(self):
        '''
        NEST produces one population spike file per virtual process.
        This function returns the corrected GIDs,
        i.e. first corrected GID and population size of each population
        '''
        fname = os.path.join(self.output_proc_prefix, 'population_GIDs.json')
        if RANK == 0:
            #load json file if it exist
            if os.path.isfile(fname):
                with open(fname, 'rb') as fp:
                    GIDs = json.load(fp)
            else:
                raw_first_gids = [self.GIDs[i][0]
                                  for i in range(len(self.X))]
                converted_first_gids = [int(1 + np.sum(self.N_X[:i]))
                                            for i in range(len(self.X))]
        
                GIDs = {}
                for i, X in enumerate(self.X):
                    GIDs[X] = [converted_first_gids[i], self.N_X[i]]
                
                print('Writing corrected population GIDs to file:')
                print('writing:\t{0}'.format(fname))
                with open(fname, 'wb') as fp:
                    json.dump(GIDs, fp)
                print('\n')
                
                print('\nGlobal corrected ids (X :\t[first id, N_X]):')
                for X, GID in list(GIDs.items()):
                    print('{0} :\t{1}'.format(X, GID))
                print('\n')
                
        else:
            GIDs = None
        
        return COMM.bcast(GIDs, root=0)    


    def merge_files(self, detector_label):
        '''
        NEST produces one population spike/voltage file per virtual process.
        This function gathers and combines them into one single file per population.
        
        Arguments
        ---------
        detector_label : str
        
        '''

        if RANK == 0:
            print('Writing spikes/voltages with corrected GIDs to file:')

        for pop_idx in range(len(self.X)):
            #parallelize on the population level
            if pop_idx % SIZE == RANK:
                files = glob.glob(os.path.join(self.output_raw_prefix,
                            detector_label + str(pop_idx) + '*'))
                gdf = []
                for f in files:
                    new_gdf = helpers.read_gdf(f)
                    # spike files: GID, spike time
                    # voltage files: GID, voltage probe time, voltage
                    for el in new_gdf:
                        # correct GID
                        el[0] = self.correct_any_GID(el[0])
                        # subtract transient from NEST simulation so that
                        # spike times start at 0.
                        el[1] -= self.transient
                        gdf.append(el)

                if 'voltages' in detector_label:
                    # sort for neuron ids (along the first axis)
                    gdf = sorted(gdf, key=lambda x: x[0])

                print('writing: %s' % os.path.join(self.output_proc_prefix,
                    detector_label + '%s.gdf' % self.X[pop_idx]))
                helpers.write_gdf(gdf, os.path.join(self.output_proc_prefix,
                    detector_label + '%s.gdf' % self.X[pop_idx]))
        print('\n')
                
        COMM.Barrier()

 
    def correct_any_GID(self, aGID):
        '''
        Needs self.GIDs and self.GIDs_corrected. Checks to which population a
        GID belongs and corrects it
        
        Arguments
        ---------
        aGID : int,
            index of a GID in population
        
        '''
        for i,X in enumerate(self.X):
            if (aGID >= self.GIDs[i][0]) and (aGID <= self.GIDs[i][1]):
                GID_corrected = self.GIDs_corrected[X][0] + aGID - self.GIDs[i][0]
                break
            else:
                continue
        return GID_corrected
  

    def get_positions(self):
        '''
        Get the neuron positions from file.
        '''
        # convert lists to a nicer format, i.e., [[L23E, L23I], []....]
        print('Loading positions from file')
        #do processing on RANK 0, but
        #consider parallel implementation, however procedure is fairly fast
        if RANK == 0:
            #attempt to reload file if it exist
            fname = os.path.join(self.output_proc_prefix, 'all_positions.h5')
            if os.path.isfile(fname):
                #hdf5 container
                all_positions = h5py.File(fname)
                
                #dict data container, fill in values
                all_positions_dict = {}
                for key, value in list(all_positions.items()):
                    all_positions_dict[key] = value.value
                
                all_positions.close()
            else:
                #position files
                files = glob.glob(os.path.join(self.output_raw_prefix,
                                               self.position_filename + '*.dat'))
                #iterate over files    
                for i, f in enumerate(files):
                    if i == 0:
                        positions = helpers.read_gdf(f)
                    else:
                        positions = np.r_[positions, helpers.read_gdf(f)]
    
                #sort according to GID
                positions = positions[positions[:, 0].argsort()]
    
                #correct GIDs
                for i, pos in enumerate(positions):
                    positions[i][0] = self.correct_any_GID(positions[i][0])
    
                #dict data container
                all_positions_dict = {}
    
                #hdf5 container
                all_positions = h5py.File(fname)
    
                #fill in values
                for i, X in enumerate(self.X):
                    X_pos = []
                    for j, pos in enumerate(positions):
                            # if GID belongs to population X
                            if (pos[0] >= self.GIDs_corrected[X][0]) and \
                                (pos[0] <= sum(self.GIDs_corrected[X])-1):
                                X_pos.append(pos.tolist())
                    all_positions_dict[X] = np.array(X_pos)
    
                    #dump to h5 if entry doesn't exist:
                    if X not in list(all_positions.keys()):
                        # each dataset has a fixed type -> here: GIDs have become floats
                        dset = all_positions.create_dataset(X,
                                                            data=X_pos,
                                                            compression='gzip',
                                                            compression_opts=2)
    
                all_positions.close()
            print('Positions loaded.')
        
        else:
            all_positions_dict = None
        return COMM.bcast(all_positions_dict, root=0)
            

    def _get_connections(self):
        '''
        Get the recorded connections from file.
        '''
        #do processing on RANK 0, but consider parallel implementation
        #if RANK == 0:
        print('Loading connections from file\n')
        files = glob.glob(os.path.join(self.output_raw_prefix,
                                       self.connection_filename + '*.dat'))
        if len(files) == 0:
            all_raw_connections = None
            pass
        else:
            #iterate over files    
            for i, f in enumerate(files):
                if i == 0:
                    connections = helpers.read_gdf(f)
                else:
                    connections = np.r_[connections, helpers.read_gdf(f)]
            
            #sort according to GID
            connections = connections[connections[:, 0].argsort()]
            
            # sort by the population of the source neuron and
            # convert lists to a nicer format
            all_raw_connections = {}
            for i,X in enumerate(self.X):
                all_raw_connections[X] = []
    
            for i in range(len(connections)):
                for j,X in enumerate(self.X):
                    if (connections[i][0] >= self.GIDs[j][0]) \
                    and (connections[i][0] <= self.GIDs[j][1]):
                        all_raw_connections[X] += [list(connections[i])]
                        break
                    else:
                        continue
            
            for i,X in enumerate(self.X):
                all_raw_connections[X] = np.array(all_raw_connections[X],
                                                  dtype=object)
            print('Connections loaded.')
            return all_raw_connections


    def get_connections(self):
        '''
        Corrects the GIDs of source and target neurons
        and writes them to file.
        '''
        
        #run only on rank 0
        if RANK == 0:
            #file with corrected connection data
            fname = os.path.join(self.output_proc_prefix, 'all_connections.h5')
            
            if os.path.isfile(fname):
                print('Getting connections with corrected GIDs from file.')

                #hdf5 container
                f = h5py.File(fname)
                
                #fill in dict
                conn_corr = {}
                for key, value in list(f.items()):
                    conn_corr[key] = value.value

                f.close()
            else:                
                print('Writing connections with corrected GIDs to file.')

                #get connections
                connections = self._get_connections()
                #hdf5 container
                f = h5py.File(fname)
           
                #if no connections were written during simulation,
                #close file, return None
                if connections == None:
                    f.close()
                    conn_corr = None
                    pass

                else:
                    # a normal copy result in overwritten values in connections
                    conn_corr = copy.deepcopy(connections)
                    for i, X in enumerate(self.X):
                        if not conn_corr[X].size: # needed for hdf on hambach, crashes if empty
                            continue
                        for j, conn in enumerate(conn_corr[X]):
                            conn_corr[X][j][0] = self.correct_any_GID(conn_corr[X][j][0])
                            conn_corr[X][j][1] = self.correct_any_GID(conn_corr[X][j][1])
                        
                        ##dump to h5 if entry doesn't exist:
                        #if X not in f.keys():
                        #dump to h5 file
                        dset = f.create_dataset(X,
                                                data=conn_corr[X].tolist(),
                                                compression='gzip', compression_opts=2)

                    f.close()
        
        else:
            conn_corr = None
        return COMM.bcast(conn_corr, root=0)


    def compute_position_hist(self, positions):
        '''
        assign the neurons to spatial bins

        Arguments
        ---------
        positions : np.ndarray
            columns are [neuron #, x-pos, y-pos]

        Returns
        -------
        pos_hist : np.ndarray
            histogram over neuron positions
        '''
        pos_hist = np.histogram2d(positions[:, 2].astype(float),
                                          positions[:, 1].astype(float),
                                  bins=[self.pos_bins, self.pos_bins])[0]
        return pos_hist.astype(int)


    def compute_time_binned_sptrains(self, X, gdf, time_bins, dtype=np.uint16):
        '''
        Compute the spike train histograms

        Arguments
        ---------        
        X : str
            population name
        gdf : np.ndarray
            colums are [GID, spike_time]
        time_bins : np.ndarray
            bin array for histogram
        dtype : type(int)
            any integer type that will fit data.

        Returns
        -------
        out : scipy.sparse.csr.csr_matrix
            sparse pop-size times time_bins.size spike rate histogram
        
        '''
        #need one extra bin
        dt = np.diff(time_bins)[0]
        time_bins_h = np.r_[time_bins, [time_bins[-1] + dt]]

        #spike train histogram, use scipy.sparse.lil_matrix
        #(row-based linked list sparse matrix).
        #lil_matrix supports slicing which coo_matrix do not, so we convert to
        #coo (COOrdinate format) later on before saving to file
        sptrains = sp.lil_matrix((self.rec_sizes[self.X.index(X)],
                                  time_bins.size), dtype=dtype)
        
        #if no spikes, return empty array
        if gdf.size == 0:
            return sptrains.tocsr()
        else:
            # get indices of the time bins to which each spike time belongs
            pos_t = np.digitize(gdf[:, 1].astype(np.float), time_bins_h[1:])
                    
            #create COO matrix 
            sptrains = sp.coo_matrix((np.ones(gdf.shape[0], dtype=dtype),
                                      (gdf[:,0]-self.GIDs_corrected[X][0], pos_t)),
                shape=sptrains.shape, dtype=dtype)
        

            # TODO: OUT OF BOUNDS IF ONLY ONE TC NEURON, CHECK IT!

            #we are not changing the matrix later, so convert to CSR_matrix
            #that supports faster iteration
            return sptrains.tocsr()
        

    def compute_binned_rates(self, binned_spcounts, pos_hist):
        '''
        Compute the spatially binned spike rates

        Arguments
        ---------
        binned_spcounts : np.ndarray
            spatiotemporally binned counts of spikes
        pos_hist : np.ndarray
            histogram over bin positions

        Returns
        -------
        binned_sprates : np.ndarray
            averaged spike rates over spatial bins

        '''
        #spatial spike rates, avoid division by zero, start TRANSIENT is removed
        binned_sprates = np.zeros(binned_spcounts.shape)
        binned_sprates[pos_hist != 0] = binned_spcounts[pos_hist != 0].astype(
            float) / pos_hist[pos_hist != 0]
        binned_sprates *= 1E3
        binned_sprates /= (self.t_sim - self.TRANSIENT)

        return binned_sprates


    def compute_binned_sprates(self, binned_sptrains, pos_hist, time_bins):
        '''
        compute the the time and position downsampled spike rates averaged over
        neurons in a spatial bin
        
        Arguments
        ---------
        binned_sptrains : np.ndarray
            time-binned spike data
        pos_hist : np.ndarray
            histogram over neuron positions
        time_bins : np.ndarray
            array over timebins
                
        Returns
        -------
        binned_sprates : np.ndarray
            time and position-binned spike rate
        
        '''
        ##avoid division by zero, flatten for array division
        pos_hist = pos_hist.flatten()
        pos_inds = np.where(pos_hist > 0)
        #pos_denom = pos_hist[pos_mask]
        
        pos_hist_diag = sp.diags(pos_hist.flatten(), 0)

        #2D-array for instantaneous rate, shaped as input
        binned_sprates = binned_sptrains.astype(float)

        #compute rate profile as:
        #(spike train per patch) / (neurons per patch) * 1000 ms/s / dt
        #hence we do not divide by duration but account for dt
        binned_sprates *= 1E3
        binned_sprates /= self.BINSIZE_TIME

        binned_sprates = binned_sprates.tolil()
        binned_sprates[pos_inds] = (binned_sprates[pos_inds].toarray().T /  pos_hist[pos_inds]).T
    
        return binned_sprates.tocsr()
    

    def compute_pos_binned_sptrains(self, positions, sptrains, dtype=np.uint8):
        '''
        compute the position-binned spike trains

        Arguments
        ---------
        positions : np.ndarray
            neuron positions
        sptrains : np.ndarray
            neuron spike trains
        dtype : type(int)
            integer type that can fit data
        
        Returns
        -------
        pos_binned_sptrains : np.ndarray
            pos x pos x spike rate 3D spike rate array
        
        
        TODO: Look into using sparse ndarray type outside of scipy.sparse
        for output

        '''
        #match position indices with spatial indices
        pos_x = np.digitize(positions[:, 1].astype(float), self.pos_bins[1:], right=True)
        pos_y = np.digitize(positions[:, 2].astype(float), self.pos_bins[1:], right=True)
        
        #2D sparse array with spatial bins flattened to 1D
        map_y, map_x = np.mgrid[0:self.pos_bins.size-1, 0:self.pos_bins.size-1]
        map_y = map_y.ravel()
        map_x = map_x.ravel()
        
        sptrains_coo = sptrains.tocoo().astype(dtype)
        nspikes = np.asarray(sptrains_coo.sum(axis=1)).flatten().astype(int)
        data = sptrains_coo.data
        col = sptrains_coo.col
        row = np.zeros(sptrains_coo.nnz)
        j = 0
        for i, n in enumerate(nspikes):
            [ind] = np.where((map_x == pos_x[i]) & (map_y==pos_y[i]))[0]
            row[j:j+n] = ind
            j += n
        
        binned_sptrains_coo = sp.coo_matrix((data, (row, col)), shape=(map_x.size, sptrains.shape[1]))
        
        try:
            assert(sptrains.sum() == binned_sptrains_coo.tocsr().sum())
        except AssertionError as ae:
            raise ae('sptrains.sum()={0} != binned_sptrains_coo.sum()={1}'.format(sptrains.sum(),
                                                                                   binned_sptrains_coo.tocsr().sum()))
        
        return binned_sptrains_coo.tocsr()
        

    def compute_pos_binned_voltages(self, X, positions, voltages, binned_positions):
        '''
        Perform a spatial binning on membrane voltage data. Compute the mean
        membrane voltage in each bin at each time step..
        '''

        # at how many time points voltages were probed
        num_time_steps = len(np.unique(voltages[:,1]))

        # one list for the voltage trace of each neuron
        voltages = np.reshape(voltages, (-1, num_time_steps, 3))

        #match position indices with spatial indices
        pos_x = np.digitize(positions[:, 1].astype(float), self.pos_bins[1:], right=True)
        pos_y = np.digitize(positions[:, 2].astype(float), self.pos_bins[1:], right=True)

        # create empty 3D array
        binned_voltages = np.zeros((self.pos_bins.size-1,
                                    self.pos_bins.size-1,
                                    num_time_steps), dtype=float)

        for i in range(len(voltages)):
            # index is neuron id - first id of this population
            idx = voltages[i][0][0] - self.GIDs_corrected[X][0]
            # fill voltage traces into array, computing the mean voltage
            # if more than one neuron is inside a bin using the position histogram
            binned_voltages[pos_x[idx], pos_y[idx], ] += voltages[i][:,2] \
                / float(binned_positions[pos_x[idx], pos_y[idx]])

        return binned_voltages


    def get_pos_sorting_array(self, positions):
        '''
        get sorting arrays that would sort GIDs according to which axis

        Arguments
        ---------
        positions : np.ndarray
            neuron positions

        Returns
        -------
        argsort : np.ndarray
            sorting array

        '''
        mssg = "sorting axis must be 'x', 'y' or None!"
        assert self.sorting_axis in ['x', 'y', None], mssg
        if self.sorting_axis == 'x':
            argsort = np.argsort(positions[:, 1])
        elif params.sorting_axis == 'y':
            argsort = np.argsort(positions[:, 2])
        else:
            argsort = np.arange(positions.shape[0])

        return argsort


    def get_neuron_pair_distances(self, positions, max_num_pairs=100):
        '''
        compute the distances between pairs of neurons within the same population.
        for a population size of N, there are N(N-1)/2 pairs,
        but only atmost 'max_num_pairs' ones are considered

        the implemented method (drawing pairs of indeces and accept them if they were
        not drawn before) is appropriate for drawing a few number of pairs randomly
        from a huge number of possible pairs.

        Arguments
        ---------
        positions: numpy.ndarray
            array with lists (id, x-position, y-position)
        max_num_pairs: int
            maximum number of pairs to be evaluated for this population
                       
        Returns
        -------
        pair_dist : numpy.ndarray
            arrays with lists (id1, id2, distance)
        '''

        pair_dist_pop = []

        # for getting all possible pairs
        #        for i in np.arange(len(positions)-1):
        #            for j in np.arange(i+1, len(positions)): # don't count pairs twice
        #                if counter < limit: # less than the maximum number of pairs
        #                    diff = positions[j][1:] - positions[i][1:]
        #                    distance = math.hypot(diff[0], diff[1]) # a bit faster than np.sqrt()...
        #                    id1 = positions[i][0]
        #                    id2 = positions[j][0]
        #                    pair_dist_pop.append([id1, id2, distance])
        #                else: 
        #                    continue

        len_pos = len(positions)
        possible_max_num_pairs = int((len_pos*(len_pos-1.))/2.)
        if max_num_pairs > possible_max_num_pairs:
            max_num_pairs = possible_max_num_pairs
        
        # maximum number of attempts to draw test pairs
        max_counter = max_num_pairs * 2

        pair_inds = [] # indeces of pairs
        uset = set()
        counter = 0 # how many test pair indeces are generated
        while len(pair_inds) < max_num_pairs and counter < max_counter:
            counter += 1
            test_pair_ind = np.sort([np.random.randint(0,len_pos-1),
                                     np.random.randint(0,len_pos-1)])
            str_pair = test_pair_ind.tostring()
            if test_pair_ind[0] == test_pair_ind[1]:
                continue
            if str_pair in uset: # if pair was not not drawn before
                continue
            else:
                pair_inds.append([test_pair_ind])
                uset.add(str_pair)

                i = test_pair_ind[0]
                j = test_pair_ind[1]
                diff = np.abs(positions[j][1:] - positions[i][1:])
                if self.pbc:
                    diff = [d%(self.extent_length/2.) for d in diff]
                distance = math.hypot(diff[0], diff[1]) # a bit faster than np.sqrt()...
                id1 = positions[i][0]
                id2 = positions[j][0]
                pair_dist_pop.append([id1, id2, distance])

        pair_dist = copy.deepcopy(pair_dist_pop)
        pair_dist.sort() # sort by first argument = first id
        
        return np.array(pair_dist)
    

    def print_firing_rates(self):
        '''
        compute and print the firing rates from spike detector output
        
        Returns
        -------
        output : list,
            each list element a length 2 list of population name and rate
        
        '''
        if RANK == 0:
            rates = []
            for i, X in enumerate(self.X):
                if self.spike_detector_output[X].size > 0:
                    rates += [float((self.spike_detector_output[X][:, 1] >= self.TRANSIENT).sum()) \
                              / self.rec_sizes[i] / (self.t_sim-self.TRANSIENT) * 1E3]
                else:
                    rates += [0]
    
            outp = list(zip(self.X, rates))
                
            print('Firing rates:')
            for X, r in outp:
                print('{0}:\t{1}'.format(X, r))
            print('\n')
            
        else:
            outp = None
        return COMM.bcast(outp, root=0)
