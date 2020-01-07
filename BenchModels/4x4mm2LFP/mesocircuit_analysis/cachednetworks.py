#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Implementation of hybridLFPy.CachedNetwork class with network populations
created using the NEST topology library.

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

import os
from glob import glob
import numpy as np
import h5py
from hybridLFPy import CachedNetwork, GDF
from mpi4py import MPI


#################################################
# Initialization of MPI stuff                   #
#################################################
COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()


#################################################
# Define derived classes for network spikes     #
#################################################

class CachedTopoNetwork(CachedNetwork):
    def __init__(self,
                 autocollect=True,
                 label_positions='brunel-py-pos',
                 **kwargs):
        '''
        
        Parameters
        ----------
        autocollect : bool
            whether or not to automatically gather gdf file output
        label_positions : str
            file prefix of position txt files
        **kwargs :
            parameters for parent class hybridLFPy.CachedNetwork
        '''
        #initialize parent class
        CachedNetwork.__init__(self, autocollect=autocollect, **kwargs)

        #set class attributes
        self.label_positions = label_positions
        
        #load positions and set them as attributes
        self.positions = {}
        for X in self.X:
            fname = os.path.join(self.spike_output_path, 'all_positions.h5')
            if os.path.isfile(fname):
                f = h5py.File(fname)
                #set positions, units from mm to mum !!!!!!!!!!!!!!!!!!!!!!!!!!!
                self.positions[X] = f[X].value[:, 1:]*1E3
                f.close()
            else:
                fnames = glob(os.path.join(self.spike_output_path,
                                          label_positions+'*{0}*txt'.format(X)))
                for i, fname in enumerate(fnames):
                    if i == 0:
                        tmp_pos = np.loadtxt(fname, dtype=object)
                    else:
                        tmp_pos = np.vstack((tmp_pos,
                                             np.loadtxt(fname, dtype=object)))
                #sorting array
                argsort = np.argsort(tmp_pos[:, 0].astype(int))                    
                
                #set positions
                self.positions[X] = tmp_pos[argsort, 1:].astype(float)

    
    def plot_raster(self, ax, xlim, x, y, pop_names=False,
                    markersize=20., alpha=1., legend=True, ):
        """
        Plot network raster plot in subplot object.
        
        
        Parameters
        ----------
        ax : `matplotlib.axes.AxesSubplot` object
            plot axes
        xlim : list
            List of floats. Spike time interval, e.g., [0., 1000.].
        x : dict
            Key-value entries are population name and neuron spike times.
        y : dict
            Key-value entries are population name and neuron gid number.
        pop_names: bool
            If True, show population names on yaxis instead of gid number.
        markersize : float
            raster plot marker size
        alpha : float in [0, 1]
            transparency of marker
        legend : bool
            Switch on axes legends.


        Returns
        -------
        None
        
        """
        for i, X in enumerate(self.X):
            ax.plot(x[X], y[X], 'o',
                markersize=markersize,
                markerfacecolor=self.colors[i],
                markeredgecolor=self.colors[i],
                alpha=alpha,
                label=X, rasterized=True,
                clip_on=True)
        
        ax.axis([xlim[0], xlim[1], 0, self.N_X.sum()])
        ax.set_ylim(ax.get_ylim()[::-1])
        ax.set_ylabel('cell id', labelpad=0)
        ax.set_xlabel('$t$ (ms)', labelpad=0)
        if legend:
            ax.legend()
        if pop_names:
            yticks = []
            yticklabels = []
            for i, X in enumerate(self.X):
                if y[X] != []:
                    yticks.append(y[X].mean())
                    yticklabels.append(self.X[i])
            ax.set_yticks(yticks)
            ax.set_yticklabels(yticklabels)

        # Add some horizontal lines separating the populations
        for X in self.X:
            if y[X].size > 0:
                ax.plot([xlim[0], xlim[1]], [y[X].max(), y[X].max()],
                    'k', lw=0.25)        
    


class CachedFixedSpikesTopoNetwork(CachedTopoNetwork):
    """
    Subclass of CachedTopoNetwork.
    
    Fake nest output, where each cell in a subpopulation spike
    simultaneously, and each subpopulation is activated at times given in
    kwarg activationtimes.
    
    
    Parameters
    ----------
    activationtimes : list of floats
        Each entry set spike times of all cells in each population
    autocollect : bool
        whether or not to automatically gather gdf file output
    **kwargs : see parent class `hybridLFPy.cachednetworks.CachedNetwork`
    
    
    Returns
    -------
    `hybridLFPy.cachednetworks.CachedFixedSpikesNetwork` object


    See also
    --------
    hybridLFPy.CachedNetwork, hybridLFPy.CachedNoiseNetwork, 

    """
    def __init__(self,
                 activationtimes=[200, 300, 400, 500, 600, 700, 800, 900, 1000],
                 filelabel='population_spikes',
                 mask=[-400, 0, -400, 0],
                 autocollect=False,
                 **kwargs):
        """
        Subclass of CachedNetwork
        
        Fake nest output, where each cell in a subpopulation spike
        simultaneously, and each subpopulation is activated at times given in
        kwarg activationtimes.
        
        Parameters
        ----------
        activationtimes : list
            Each entry set spike times of all cells in each population
        filelabel : str
            file prefix for spike time files
        mask : list
            Coordinates [xmin, xmax, ymin, ymax]. Only neurons within mask will
            be assigned a fixed spike time. 
        autocollect : bool
            whether or not to automatically gather gdf file output
        **kwargs : see parent class `hybridLFPy.cachednetworks.CachedNetwork`
        
        
        Returns
        -------
        `hybridLFPy.cachednetworks.CachedFixedSpikesNetwork` object


        See also
        --------
        CachedNetwork, CachedNoiseNetwork, 
        
        """
        
        CachedTopoNetwork.__init__(self, autocollect=autocollect, **kwargs)

        # Set some attributes
        self.activationtimes = activationtimes
        self.filelabel = filelabel
        self.mask = mask
        
        if len(activationtimes) != len(self.N_X):
            raise Exception('len(activationtimes != len(N_X))')

        """ Create a dictionary of nodes with proper layernames
         self.nodes = {}.
        """
        self.mask_nodes = {}

        if RANK == 0:
            for i, N in enumerate(self.N_X):
                nodes = self.nodes[self.X[i]]
                
                #apply mask
                pos = self.positions[self.X[i]]
                inds = (pos[:, 0] >= self.mask[0]) & (pos[:, 0] <= self.mask[1]) & (pos[:, 1] >= self.mask[2]) & (pos[:, 1] <= self.mask[3])
                self.mask_nodes[self.X[i]] = nodes[inds]
                
                cell_spt = list(zip(nodes[inds], [self.activationtimes[i]
                                  for x in range(nodes[inds].size)]))
                cell_spt = np.array(cell_spt, dtype=[('a', int), ('b', float)])

                np.savetxt(os.path.join(self.spike_output_path,
                                        self.filelabel + '_{0}.gdf'.format(self.X[i])),
                           cell_spt, fmt=['%i', '%.1f'])

        # Resync
        COMM.barrier()
        

        # Collect the gdf files
        self.collect_gdf()


    def collect_gdf(self):
        """
        Collect the gdf-files from network sim in folder `spike_output_path`
        into sqlite database, using the GDF-class.
        
        
        Parameters
        ----------
        None
        
        
        Returns
        -------
        None
        
        """
        # Resync
        COMM.Barrier()

        # Raise Exception if there are no gdf files to be read
        if len(glob(os.path.join(self.spike_output_path,
                                      self.filelabel + '*.'+ self.ext))) == 0:
            raise Exception('path to files contain no gdf-files!')

        #create in-memory databases of spikes
        if not hasattr(self, 'dbs'):
            self.dbs = {}
        
        for X in self.X:
            db = GDF(os.path.join(self.dbname),
                     debug=True, new_db=True)
            db.create(re=os.path.join(self.spike_output_path,
                                      '{0}*{1}*{2}'.format(self.filelabel, X,
                                                           self.ext)),
                      index=True)
            self.dbs.update({
                    X : db
                })
      
        COMM.Barrier()


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    
