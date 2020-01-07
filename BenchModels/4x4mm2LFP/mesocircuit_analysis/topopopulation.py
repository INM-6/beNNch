#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Implementation of hybrid scheme population class with distance-dependent
connectivity.

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
import numpy as np
from time import time 
from mpi4py import MPI
from hybridLFPy import Population
from .helperfun import _getSpCell, _calc_radial_dist_to_cell, _fetchSpCells, _get_all_SpCells


################# Initialization of MPI stuff ##################################
COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()


class TopoPopulation(Population):
    def __init__(self,
                 topology_connections = {
                    'EX': {
                        'edge_wrap': True,
                        'extent': [4000., 4000.],
                        'allow_autapses' : True,
                        'kernel': {'exponential': {
                            'a': 1., 'c': 0.0, 'tau': 300.}},
                        'mask': {'circular': {'radius': 2000.}},
                        'delays' : {
                            'linear' : {
                                'c' : 1.,
                                'a' : 2.
                                }
                            },
                        },
                     'IN': {
                        'edge_wrap': True,
                        'extent': [4000., 4000.],
                        'allow_autapses' : True,
                        'kernel': {'exponential': {
                            'a': 1., 'c': 0.0, 'tau': 300.}},
                        'mask': {'circular': {'radius': 2000.}},
                        'delays' : {
                            'linear' : {
                                'c' : 1.,
                                'a' : 2.
                                }
                            },                            
                        },
                    },
                 **kwargs):
        '''
        Initialization of class TopoPopulation, for dealing with networks
        created using the NEST topology library (distance dependent
        connectivity). 
        
        Inherited of class hybridLFPy.Population
        
        Arguments
        ---------
        topology_connections : dict
            nested dictionary with topology-connection parameters for each
            presynaptic population
        
        
        Returns
        -------
        object : populations.TopoPopulation
            population object with connections, delays, positions, simulation
            methods
            
        See also
        --------
        hybridLFPy.Population
        
        '''
        #set networkSim attribute so that monkey-patched methods can work
        self.networkSim = kwargs['networkSim']
        self.topology_connections = topology_connections
        
        #initialize parent class
        Population.__init__(self, **kwargs)
        
        #set class attributes


    def get_all_synDelays(self):
        """
        Create and load arrays of connection delays per connection on this rank
        
        Get random normally distributed synaptic delays,
        returns dict of nested list of same shape as SpCells.
    
        Delays are rounded to dt.
        
        This function takes no kwargs.

        
        Parameters
        ----------
        None
        

        Returns
        -------
        dict
            output[cellindex][populationname][layerindex]`, np.array of
            delays per connection.
        
        
        See also
        --------
        numpy.random.normal
        
        """
        tic = time() #timing

        #ok then, we will draw random numbers across ranks, which have to
        #be unique per cell. Now, we simply record the random state,
        #change the seed per cell, and put the original state back below.
        randomstate = np.random.get_state()

        #container
        delays = {}

        for cellindex in self.RANK_CELLINDICES:
            #set the random seed on for each cellindex
            np.random.seed(self.POPULATIONSEED + cellindex + 2*self.POPULATION_SIZE)

            delays[cellindex] = {}
            for j, X in enumerate(self.X):
                delays[cellindex][X] = []                    
                if 'delays' not in list(self.topology_connections[X][self.y].keys()):
                    #old behaviour, draw delays from normal distribution
                    for i, k in enumerate(self.k_yXL[:, j]):
                        loc = self.synDelayLoc[j]
                        loc /= self.dt
                        scale = self.synDelayScale[j]
                        if scale is not None:
                            scale /= self.dt
                            delay = np.random.normal(loc, scale, k).astype(int)
                            inds = delay < 1
                            while np.any(inds):
                                delay[inds] = np.random.normal(loc, scale,
                                                            inds.sum()).astype(int)
                                inds = delay < 1
                            delay = delay.astype(float)
                            delay *= self.dt
                        else:
                            delay = np.zeros(k) + self.synDelayLoc[j]
                        delays[cellindex][X] += [delay]
                else:
                    topo_conn = self.topology_connections[X][self.y]['delays']
                    if 'linear' in list(topo_conn.keys()):
                        #ok, we're using linear delays,
                        #delay(r) = a*r + c
                        a = topo_conn['linear']['a']
                        c = topo_conn['linear']['c']
                        
                        #radial distance to all cells
                        r = _calc_radial_dist_to_cell(self.pop_soma_pos[cellindex]['x'],
                                                          self.pop_soma_pos[cellindex]['y'],
                                                          self.networkSim.positions[X],
                                                          self.topology_connections[X][self.y]['extent'][0],
                                                          self.topology_connections[X][self.y]['extent'][1],
                                                          self.topology_connections[X][self.y]['edge_wrap'])
                        
                        #get presynaptic unit GIDs for connections
                        i0 = self.networkSim.nodes[X][0]
                        for i, k in enumerate(self.k_yXL[:, j]):
                            x = self.SpCells[cellindex][X][i]
                            if self.synDelayScale[j] is not None:
                                scale = self.synDelayScale[j]
                                delay = np.random.normal(0, scale, k)
                                #avoid zero or lower delays
                                inds = delay < self.dt-c
                                while np.any(inds):
                                    #print inds.size
                                    delay[inds] = np.random.normal(0, scale,
                                                                inds.sum())
                                    inds = delay < self.dt-c
                            else:
                                delay = np.zeros(x.size)
                            
                            #add linear dependency
                            delay += r[x-i0]*a + c

                            #round delays to nearest dt
                            delay /= self.dt
                            delay = np.round(delay) * self.dt

                            #fill in values
                            delays[cellindex][X] += [delay]
                    else:
                        raise NotImplementedError('{0} delay not implemented'.format(list(topo_conn.keys())[0]))

        #reset the random number generator
        np.random.set_state(randomstate)

        if RANK == 0:
            print(('found delays in %.2f seconds' % (time()-tic)))

        return delays
            
        
    def get_all_SpCells(self):
        """
        For each postsynaptic cell existing on this RANK, load or compute
        the presynaptic cell index for each synaptic connection to the
        postsynaptic cell, using distance-dependent connectivity

        This function takes no kwargs.

        
        Parameters
        ----------
        None
        
        
        Returns
        -------
        SpCells : dict
            `output[cellindex][populationname][layerindex]`, np.array of
            presynaptic cell indices.
        
        
        See also
        --------
        Population.fetchSpCells, TopoPopulation.fetchSpCells
        
        """
        tic = time() #timing
        
        
        #ok then, we will draw random numbers across ranks, which have to
        #be unique per cell. Now, we simply record the random state,
        #change the seed per cell, and put the original state back below.
        randomstate = np.random.get_state()
        
        #set the random seed on for first cellindex on RANK
        if self.RANK_CELLINDICES.size > 0:
            np.random.seed(self.POPULATIONSEED + self.RANK_CELLINDICES[0] + self.POPULATION_SIZE)
        else:
            pass #no cells should be drawn here anyway.



        
        SpCells = _get_all_SpCells(self.RANK_CELLINDICES,
                                   self.X,
                                   self.y,
                                   self.pop_soma_pos,
                                   self.networkSim.positions,
                                   self.topology_connections,
                                   self.networkSim.nodes,
                                   self.k_yXL)

        #reset the random number generator
        np.random.set_state(randomstate)

        if RANK == 0:
            print(('found presynaptic cells in %.2f seconds' % (time()-tic)))

        #resync
        COMM.Barrier()
    
        return SpCells


    def set_pop_soma_pos(self):
        """
        Set `pop_soma_pos` using draw_rand_pos().
        
        This method takes no keyword arguments.
        
        
        Parameters
        ----------
        None
        
        
        Returns
        -------
        np.ndarray
            (x,y,z) coordinates of each neuron in the population
        
        
        See also
        --------
        TopoPopulation.draw_rand_pos
        
        """
        if RANK == 0:
            pop_soma_pos = self.draw_rand_pos(**self.populationParams)
        else:
            pop_soma_pos = None
        return COMM.bcast(pop_soma_pos, root=0)
    

    def draw_rand_pos(self, z_min, z_max, position_index_in_Y, **args):
        """
        Draw cell positions from the CachedNetwork object (x,y), but with
        randomized z-position between z_min and z_max.
        Returned argument is a list of dicts
        [{'x' : float, 'y' : float, 'z' : float}, {...}].
        
        
        Parameters
        ----------        
        z_min : float
            Lower z-boundary of population.
        z_max : float
            Upper z-boundary of population.
        position_index_in_Y : list
            parent population Y of cell type y in Y, first index for slicing pos
        **args : keyword arguments
            Additional inputs that is being ignored.


        Returns
        -------
        soma_pos : list
            List of dicts of len population size
            where dict have keys x, y, z specifying
            xyz-coordinates of cell at list entry `i`.
        
        
        See also
        --------
        TopoPopulation.set_pop_soma_pos
           
        """

        print("assess somatic locations: ")
        Y, ind = position_index_in_Y
        z = np.random.rand(self.POPULATION_SIZE)*(z_max - z_min) + z_min
        soma_pos = [{
            'x' : x,
            'y' : y,
            'z' : z[k]
            } for k, (x, y) in enumerate(self.networkSim.positions[Y][ind:ind+z.size])]

        print('done!')

        return soma_pos


if __name__ == '__main__':
    import doctest
    doctest.testmod()
