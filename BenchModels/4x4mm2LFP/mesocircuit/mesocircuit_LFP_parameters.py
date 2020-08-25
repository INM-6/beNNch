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
import sys
try: # prepending to PYTHONPATH doesn't work on JURECA, for example matplotlib must be updated to version 2* to make correct plots
    sys.path.remove('/usr/local/software/jureca/Stages/2016a/software/SciPy-Stack/2016a-intel-para-2016a-Python-2.7.11/lib/python2.7/site-packages/matplotlib-1.5.1-py2.7-linux-x86_64.egg')
except ValueError:
    pass 
import os
if 'DISPLAY' not in os.environ:
    import matplotlib
    matplotlib.use('Agg')
import sys
import numpy as np
from mesocircuit_analysis import helpers
import parameterspace_control as psc
import dataset_analysis as dsa
from NeuroTools.parameters import ParameterSet
import json



####################################
# HELPER FUNCTIONS                 #
####################################

flattenlist = lambda lst: sum(sum(lst, []),[])


def kill_sli_in_nested_dict(slidict):
    '''
    rebuild level 2 nested input dictionary and remove slashes in key names
    '''
    newdict = dict()
    for key0, value0 in list(slidict.items()):
        key1 = key0.replace('/', '')
        newdict.update({key1 : {}})
        for key, val in list(value0.items()):
            newdict[key1].update({key.replace('/', '') : val})
    return newdict


####################################
# SPATIAL CONNECTIVITY EXTRACTION  #
####################################

'''
Include functions that extract and process information from binzegger.json here
'''


def get_F_y(fname='binzegger_connectivity_table.json', y=['p23']): 
    '''
    Extract frequency of occurrences of those cell types that are modeled.
    The data set contains cell types that are not modeled (TCs etc.)
    The returned percentages are renormalized onto modeled cell-types, i.e. they sum up to 1 
    '''
    # Load data from json dictionary
    f = open(fname,'r')
    data = json.load(f)
    f.close()
    
    occurr = []
    for cell_type in y:
        occurr += [data['data'][cell_type]['occurrence']]
    return list(np.array(occurr)/np.sum(occurr)) 


def get_L_yXL(fname, y, x_in_X, L):
    '''
    compute the layer specificity, defined as:
    ::
    
        L_yXL = k_yXL / k_yX
    '''
    def _get_L_yXL_per_yXL(fname, x_in_X, X_index,
                                  y, layer):
        # Load data from json dictionary
        f = open(fname, 'r')
        data = json.load(f)
        f.close()
        
        
        # Get number of synapses
        if layer in [str(key) for key in list(data['data'][y]['syn_dict'].keys())]:
            #init variables
            k_yXL = 0
            k_yX = 0
            
            for x in x_in_X[X_index]:
                p_yxL = data['data'][y]['syn_dict'][layer][x] / 100.
                k_yL = data['data'][y]['syn_dict'][layer]['number of synapses per neuron']
                k_yXL += p_yxL * k_yL
                
            for l in [str(key) for key in list(data['data'][y]['syn_dict'].keys())]:
                for x in x_in_X[X_index]:
                    p_yxL = data['data'][y]['syn_dict'][l][x] / 100.
                    k_yL = data['data'][y]['syn_dict'][l]['number of synapses per neuron']
                    k_yX +=  p_yxL * k_yL
            
            if k_yXL != 0.:
                return k_yXL / k_yX
            else:
                return 0.
        else:
            return 0.


    #init dict
    L_yXL = {}

    #iterate over postsynaptic cell types
    for y_value in y:
        #container
        data = np.zeros((len(L), len(x_in_X)))
        #iterate over lamina
        for i, Li in enumerate(L):
            #iterate over presynapse population inds
            for j in range(len(x_in_X)):
                data[i][j]= _get_L_yXL_per_yXL(fname, x_in_X,
                                                          X_index=j,
                                                          y=y_value,
                                                          layer=Li)
        L_yXL[y_value] = data

    return L_yXL


def get_T_yX(fname, y, y_in_Y, x_in_X, F_y):
    '''
    compute the cell type specificity, defined as:
    ::
    
        T_yX = K_yX / K_YX
            = F_y * k_yX / sum_y(F_y*k_yX) 
    
    
    '''
    def _get_k_yX_mul_F_y(y, y_index, X_index):
        # Load data from json dictionary
        f = open(fname, 'r')
        data = json.load(f)
        f.close()
    
        #init variables
        k_yX = 0.
        
        for l in [str(key) for key in list(data['data'][y]['syn_dict'].keys())]:
            for x in x_in_X[X_index]:
                p_yxL = data['data'][y]['syn_dict'][l][x] / 100.
                k_yL = data['data'][y]['syn_dict'][l]['number of synapses per neuron']
                k_yX +=  p_yxL * k_yL
        
        return k_yX * F_y[y_index]


    #container
    T_yX = np.zeros((len(y), len(x_in_X)))
    
    #iterate over postsynaptic cell types
    for i, y_value in enumerate(y):
        #iterate over presynapse population inds
        for j in range(len(x_in_X)):
            k_yX_mul_F_y = 0
            for k, yy in enumerate(sum(y_in_Y, [])):                
                if y_value in yy:
                    for yy_value in yy:
                        ii = np.where(np.array(y) == yy_value)[0][0]
                        k_yX_mul_F_y += _get_k_yX_mul_F_y(yy_value, ii, j)
            
            
            if k_yX_mul_F_y != 0:
                T_yX[i, j] = _get_k_yX_mul_F_y(y_value, i, j) / k_yX_mul_F_y
            
    return T_yX


################################################################################
# Parameter class method lifted from hybrid scheme paper work                  #
################################################################################

class general_params(object):
    def __init__(self, PSET, analysis):
        '''
        class collecting general model parameters
        
        Arguments
        ---------
        PSET : NeuroTools.parameters.ParameterSet
            Parameterset loaded from file
        analysis : dataset_analysis.NetworkAnalysis
            NetworkAnalysis object
        '''
        
        #set class attributes
        self.PSET = PSET
        self.analysis = analysis
    
    ####################################
    #                                  #
    #                                  #
    #        MODEL PARAMETERS          #
    #                                  #
    #                                  #
    ####################################



        ####################################
        # POPULATIONS                      #
        ####################################

        #GID info for each population
        self.GIDs = analysis.get_GIDs()

        # number of neurons in each population (unscaled)
        #self.full_scale_num_neurons = PSET.mm2_neuron_numbers
        self.full_scale_num_neurons = PSET.full_num_neurons

        # population names
        self.X = ['TC', 'L23E', 'L23I', 'L4E', 'L4I', 'L5E', 'L5I', 'L6E', 'L6I']
        self.Y = self.X[1:]

        # TC and cortical population sizes in one list ordered
        # as class attribute X
        self.N_X = np.r_[analysis.N_X[-1], analysis.N_X[:-1]]
          
      
        ####################################
        # CELL-TYPE PARAMETERS             #
        ####################################
        
        # Note that these parameters are only relevant for the point-neuron
        # network in case one wants to calculate depth-resolved cell-type
        # specific input currents

        # point to .json connectivity table file
        self.connectivity_table = os.path.join('binzegger',
                                               'binzegger_connectivity_table.json')

        #list of cell type names used in this script
        #names of every post-syn pop layer
        self.y_in_Y = [
                [['p23'],['b23','nb23']],
                [['p4','ss4(L23)','ss4(L4)'],['b4','nb4']],
                [['p5(L23)','p5(L56)'],['b5','nb5']],
                [['p6(L4)','p6(L56)'],['b6','nb6']]]

        self.y = flattenlist(self.y_in_Y)
        
        #need presynaptic cell type to population mapping
        self.x_in_X = [['TCs', 'TCn']] + sum(self.y_in_Y, [])
        
        
        #map the pre-synaptic populations to the post-syn populations
        self.mapping_Yy = list(zip(
                  ['L23E', 'L23I', 'L23I',
                   'L4E', 'L4E', 'L4E', 'L4I', 'L4I',
                   'L5E', 'L5E', 'L5I', 'L5I',
                   'L6E', 'L6E', 'L6I', 'L6I'],
                  self.y))

        # Frequency of occurrence of each cell type (F_y); 1-d array 
        self.F_y = get_F_y(fname=self.connectivity_table, y=self.y)

        # Relative frequency of occurrence of each cell type within its population (F_{y,Y})
        self.F_yY = [[get_F_y(fname=self.connectivity_table, y=y) for y in Y] for Y in self.y_in_Y]
        
        # Number of neurons of each cell type (N_y); 1-d array
        self.N_y =  np.array([self.full_scale_num_neurons[layer][pop] * self.F_yY[layer][pop][k] \
                                    for layer, array in enumerate(self.y_in_Y) \
                                        for pop, cell_types in enumerate(array) \
                                            for k, _ in enumerate(cell_types)])
        self.N_y *= PSET.density
        self.N_y = self.N_y.astype(int)

        
        # number of synapses in total, only intra-cortical connections
        # K_YX = np.array(PSET.full_num_synapses)
        # #number of synapses in total, only thalamocortical connections as in Potjans&Diesmann 2012
        # C_TC = np.array(sum(PSET.mm2_conn_probs_th, [])) #len 8 array
        # K_Y_TC = np.log(1. - C_TC) / np.log(1. - 1./(PSET.n_thal_per_mm2 *
        #                                              PSET.extent_length**2 *
        #                                              np.array(sum(self.full_scale_num_neurons, []))))
        # #concatenate number os synapses of TC and cortical populations
        # K_YX = np.c_[K_Y_TC, K_YX]

        #concatenate number of synapses of TC and cortical populations
        K_YX = np.c_[np.array(PSET.full_num_synapses_th),
                     np.array(PSET.full_num_synapses)].astype(float)
        
        
        #Scale the number of synapses according to network density parameter
        K_YX *= PSET.density
        K_YX = K_YX.astype(int)
        
        
        #spatial connection probabilites on each subpopulation
        #Each key must correspond to a subpopulation like 'L23E' used everywhere else,
        #each array maps thalamic and intracortical connections.
        #First column is thalamic connections, and the rest intracortical,
        #ordered like 'L23E', 'L23I' etc., first row is normalised probability of
        #connection withing L1, L2, etc.;
        self.L_yXL = get_L_yXL(fname = self.connectivity_table,
                                         y = self.y,
                                         x_in_X = self.x_in_X,
                                         L = ['1','23','4','5','6'])
        
        #compute the cell type specificity
        self.T_yX = get_T_yX(fname=self.connectivity_table, y=self.y,
                             y_in_Y=self.y_in_Y, x_in_X=self.x_in_X,
                             F_y=self.F_y)
        
        #assess relative distribution of synapses for a given celltype
        self.K_yXL = {}
        for i, (Y, y) in enumerate(self.mapping_Yy):
            #fill in K_yXL (layer specific connectivity)
            self.K_yXL[y] = (self.T_yX[i, ] * K_YX[np.array(self.Y)==Y, ] * self.L_yXL[y]).astype(int)
        
        #number of incoming connections per cell type per layer per cell 
        self.k_yXL = {}
        for y, N_y in zip(self.y, self.N_y):
            self.k_yXL.update({y : (self.K_yXL[y] / N_y).astype(int)})
        
        
    ####################################
    #                                  #
    #                                  #
    #        MODEL PARAMETERS          #
    #                                  #
    #                                  #
    ####################################

       
        ####################################
        # SCALING (DENSITY not volume)     #
        ####################################  

        self.SCALING = 1.0
        
  
        ####################################
        # MORPHOLOGIES                     #
        ####################################

        # list of morphology files with default location, testing = True
        # will point to simplified morphologies
        # self.PATH_m_y = 'morphologies'
        # self.m_y = [Y + '_' + y + '.hoc' for Y, y in self.mapping_Yy]
        self.PATH_m_y = 'morphologies'
        self.m_y = [Y + '_' + y + '.hoc' for Y, y in self.mapping_Yy]

        ####################################
        # CONNECTION WEIGHTS               #
        ####################################
        
        # compute the synapse weight from fundamentals of exp synapse LIF neuron
        self.J = self._compute_J()
        
        # set up matrix containing the synapse weights between any population X
        # and population Y, including exceptions for certain connections
        J_YX = np.array(self.PSET.weight_mod_facts) * self.J
        J_YX = np.c_[np.ones(8)*self.PSET.th_weight_scale*self.J, J_YX]
        
        # set up matrix of unit synapse weights for all possible connections
        self.J_unit = 1 
        J_YX_unit = np.zeros_like(J_YX) + self.J_unit
        
        # extrapolate weights between populations X and
        # cell type y in population Y
        self.J_yX = {}
        self.J_yX_unit = {}
        for Y, y in self.mapping_Yy:
            [i] = np.where(np.array(self.Y) == Y)[0]
            self.J_yX.update({y : J_YX[i, ]})
            self.J_yX_unit.update({y : J_YX_unit[i, ]})
        
    
        ####################################
        # GEOMETRY OF CORTICAL LAYERS      #
        ####################################
        
        # set the boundaries of each layer, L1->L6,
        # and mean depth of soma layers
        self.layerBoundaries = np.array([[     0.0,   -81.6],
                                          [  -81.6,  -587.1],
                                          [ -587.1,  -922.2],
                                          [ -922.2, -1170.0],
                                          [-1170.0, -1491.7]])
        
        # assess depth of each 16 subpopulation
        self.depths = self._calcDepths()
        
        # make a nice structure with data for each subpopulation
        self.y_zip_list = list(zip(self.y, self.m_y,
                            self.depths, self.N_y))



        ##############################################################
        # POPULATION PARAMS (cells, population, synapses, electrode) #
        ##############################################################


        # Global LFPy.Cell-parameters, by default shared between populations
        # Some passive parameters will not be fully consistent with LIF params
        self.cellParams = {
            'cm' : 1.0,
            'Ra' : 150,
            'v_init' : PSET.model_params['/E_L'],
            'passive' : True,
            'passive_parameters' : {'g_pas' : 1. / (PSET.model_params['/tau_m'] * 1E3), #assume cm=1
                                    'e_pas' : PSET.model_params['/E_L']},
            'nsegs_method' : 'lambda_f',
            'lambda_f' : 100,
            'dt' : PSET.dt,
            'tstart' : 0,
            'tstop' : PSET.t_sim,
            'verbose' : False,
        }
        

        # layer specific LFPy.Cell-parameters as nested dictionary
        self.yCellParams = self._yCellParams()
        
        
        # set the axis of which each cell type y is randomly rotated,
        # SS types and INs are rotated around both x- and z-axis
        # in the population class, while P-types are
        # only rotated around the z-axis
        self.rand_rot_axis = {}
        for y, _, _, _ in self.y_zip_list:
            #identify pyramidal cell populations:
            if y.rfind('p') >= 0:
                self.rand_rot_axis.update({y : ['z']})
            else:
                self.rand_rot_axis.update({y : ['x', 'z']})
        
        
        # additional simulation kwargs, see LFPy.Cell.simulate() docstring
        self.simulationParams = {}
        
                
        # a dict setting the number of cells N_y and geometry
        # of cell type population y
        self.populationParams = {}
        for y, _, depth, N_y in self.y_zip_list:
            self.populationParams.update({
                y : {
                    'number' : int(N_y*self.SCALING),
                    'z_min' : depth - 25,
                    'z_max' : depth + 25,
                }
            })

        # explicitly set the first neuron position index for each cell type y by
        # distributing from the proper network neuron positions
        Y0 = None
        count = 0
        for i, (Y, y) in enumerate(self.mapping_Yy):
            if Y != Y0:
                count = 0
            self.populationParams[y].update(dict(position_index_in_Y=[Y,
                                                                      count]))
            count += self.N_y[i]
            Y0 = Y
            

        # Set up cell type specific synapse parameters in terms of synapse model
        # and synapse locations
        self.synParams = {}
        for y in self.y:
            if y.rfind('p') >= 0:
                #pyramidal types have apical dendrite sections
                section = ['apic', 'dend']
            else:
                #other cell types do not
                section = ['dend']

            self.synParams.update({
                y : {
                    'syntype' : 'ExpSynI',  #current based exponential synapse
                    'section' : section,
                    # 'tau' : PSET.model_params["/tau_syn_ex"],
                },
            })


        # set up dictionary of synapse time constants specific to each
        # postsynaptic cell type and presynaptic population, such that
        # excitatory and inhibitory time constants can be varied. 
        self.tau_yX = {}
        for y in self.y:
            self.tau_yX.update({
                y : [PSET.model_params["/tau_syn_in"] if 'I' in X else
                     PSET.model_params["/tau_syn_ex"] for X in self.X]
            })


        #set parameters for topology connections with spatial parameters
        #converted to units of mum from mm.
        conn_profiles = np.c_[np.array(PSET.conn_profiles_th).flatten(),
                              PSET.conn_profiles]
        decay_params = np.c_[np.array(PSET.decay_params_th).flatten(),
                             PSET.decay_params]
        mask_radii = np.c_[np.array(PSET.mask_radii_th).flatten(),
                           PSET.mask_radii]
        
        
        #container dictionary
        self.topology_connections = {}
        for i, X in enumerate(self.X):
            self.topology_connections[X] = {}
            if 'I' in X:
                a = PSET.prop_speeds[1]*1E-3    #ms/mm -> ms/mum conversion
                c = PSET.min_delays[1]
            else:
                a = PSET.prop_speeds[0]*1E-3    #ms/mm -> ms/mum conversion
                c = PSET.min_delays[0]
                
            for Y, y in self.mapping_Yy:
                [j] = np.where(np.array(self.Y) == Y)[0]
                
                self.topology_connections[X][y] = dict(
                    extent = [PSET.extent_length*1E3,
                              PSET.extent_length*1E3],
                    edge_wrap = PSET.pbc,
                    allow_autapses = True,
                    mask = dict(
                        circular = dict(radius = mask_radii[j, i]*1E3)
                    ),
                    delays = dict(linear = dict(c = c,
                                                a = a)
                    )
                )
                #we may have different kernel types for different connections
                if conn_profiles[j, i] == 'gaussian':
                    self.topology_connections[X][y].update({
                        'kernel' : {
                            conn_profiles[j, i] : dict(
                                p_center = 1.,
                                mean = 0.,
                                c = 0.,
                                sigma = decay_params[j, i]*1E3 #mm -> mum unit conversion
                                )
                            }
                        })
                elif conn_profiles[j, i] == 'exponential':
                    self.topology_connections[X][y].update({
                        'kernel' : {
                            conn_profiles[j, i] : dict(
                                a = 1,
                                c = 0,
                                tau = decay_params[j, i]*1E3 #mm -> mum unit conversion
                                )
                            }
                        })
                elif conn_profiles[j, i] == 'random':
                    self.topology_connections[X][y].update({
                        'kernel' : conn_profiles[j, i]
                        })
                else:
                    raise NotImplementedError, 'connection profile %s not implemented' % conn_profiles[j, i]


    ############################################################################
    # CLASS METHODS PERFORMING CERTAIN CALCULATIONS                            #
    ############################################################################

    def _compute_J(self):
        '''
        Compute the current amplitude corresponding to the exponential
        synapse model PSP amplitude
        
        Derivation using sympy:
        ::
            from sympy import *
            #define symbols
            t, tm, Cm, ts, Is, Vmax = symbols('t tm Cm ts Is Vmax')
            
            #assume zero delay, t >= 0
            #using eq. 8.10 in Sterrat et al
            V = tm*ts*Is*(exp(-t/tm) - exp(-t/ts)) / (tm-ts) / Cm
            print 'V = %s' % V
            
            #find time of V == Vmax
            dVdt = diff(V, t)
            print 'dVdt = %s' % dVdt
            
            [t] = solve(dVdt, t)
            print 't(t@dVdT==Vmax) = %s' % t
            
            #solve for Is at time of maxima
            V = tm*ts*Is*(exp(-t/tm) - exp(-t/ts)) / (tm-ts) / Cm
            print 'V(%s) = %s' % (t, V)
            
            [Is] = solve(V-Vmax, Is)
            print 'Is = %s' % Is
        
        resulting in:
        ::
            Cm*Vmax*(-tm + ts)/(tm*ts*(exp(tm*log(ts/tm)/(tm - ts))
                                     - exp(ts*log(ts/tm)/(tm - ts))))
        
        '''
        #LIF params
        tm = self.PSET.model_params['/tau_m']
        Cm = self.PSET.model_params['/C_m']
        
        #synapse
        ts = self.PSET.model_params['/tau_syn_ex']
        Vmax = self.PSET.PSP_e
        
        #max current amplitude
        J = Cm*Vmax*(-tm + ts)/(tm*ts*(np.exp(tm*np.log(ts/tm)/(tm - ts))
                                     - np.exp(ts*np.log(ts/tm)/(tm - ts))))
        
        #unit conversion pF*mV -> nA
        J *= 1E-3
        
        return J


    def _calcDepths(self):
        '''
        return the cortical depth of each subpopulation
        '''
        depths = self.layerBoundaries.mean(axis=1)[1:]

        depth_y = []
        for y in self.y:
            if y in ['p23', 'b23', 'nb23']:
                depth_y = np.r_[depth_y, depths[0]]
            elif y in ['p4', 'ss4(L23)', 'ss4(L4)', 'b4', 'nb4']:
                depth_y = np.r_[depth_y, depths[1]]
            elif y in ['p5(L23)', 'p5(L56)', 'b5', 'nb5']:
                depth_y = np.r_[depth_y, depths[2]]
            elif y in ['p6(L4)', 'p6(L56)', 'b6', 'nb6']:
                depth_y = np.r_[depth_y, depths[3]]
            else:
                raise Exception('this aint right')
                
        return depth_y


    def _yCellParams(self):
        '''
        Return dict with parameters for each population.
        The main operation is filling in cell type specific morphology
        '''
        #cell type specific parameters going into LFPy.Cell        
        yCellParams = {}
        for layer, morpho, _, _ in self.y_zip_list:
            yCellParams.update({layer : self.cellParams.copy()})
            yCellParams[layer].update({
                'morphology' : os.path.join(self.PATH_m_y, morpho),
            })
        return yCellParams


################################################################################
## HYBRID SCHEME PARAMETERS
################################################################################

def get_LFP_parameters(parameter_set_file, ps_id):
    '''
    get parameter set for LFP sim corresponding to parameter set file
    
    Arguments
    ---------
    parameter_set_file : path
        full path to parameter set file
    ps_id : str
        unique identifier

    Returns
    -------
    PS : ParameterSpace
        neurotools.parameters.ParameterSpace object
    PSET : ParameterSet
        neurotools.parameters.ParameterSet object
    '''        
    #get corresponding parameterSet dictionary
    #iterate over different parameterspaces in dict
    PSET = None
    for key, PS in list(psc.ParameterSpaces.items()):
        for paramset in PS.iter_inner():
            #unique id for each parameter set, constructed from the parameset dict
            #converted to a sorted list of tuples
            id = helpers.get_unique_id(paramset)
            if ps_id == id:
                PSET = ParameterSet(os.path.join(psc.parameterset_dest, 
                                                  '{0}.pset'.format(ps_id)))
                #PSET = ParameterSet(paramset)
    
    if PSET is None:
        try:
            Warning('PSET is None, attempt loading .pset-file')
            PSET = ParameterSet(os.path.join(psc.parameterset_dest, 
                                                  '{0}.pset'.format(ps_id)))
        except NameError as ne:
            raise ne('parameterset id {} not in parameterspace.py or {}.pset file is not existing'.format(ps_id, ps_id))
        

    #set up the analysis class instance
    analysis = dsa.get_network_analysis_object(parameter_set_file, ps_id)
    
    #Set up class object with layer specificity of connections from
    #hybrid scheme implementation.
    
    #It is a convoluted way of doing it, but calculations of layer specificity of
    #connectsion and such rely on a bunch of intermediate calculations, and it
    #is not nice to do so when creating the ParameterSet object below. We'll use
    #it here as a placeholder of different parameters
    PARAMS = general_params(PSET, analysis)
    
    
    ################################################################################
    #Set up main parameters for hybrid scheme methods using the                    #
    #NeuroTools.parameters.ParameterSet class                                      #
    ################################################################################
    
    
    #set up file destinations differentiating between certain output
    PS = ParameterSet(dict(
        #Main folder of simulation output
        savefolder = os.path.join(psc.output_proc_prefix, ps_id),
        
        #make a local copy of main files used in simulations
        sim_scripts_path = os.path.join(psc.output_proc_prefix, ps_id,
                                        'sim_scripts'),
        
        #destination of single-cell output during simulation
        cells_path = os.path.join(psc.output_proc_prefix, ps_id, 'cells'),
        
        #destination of cell- and population-specific signals, i.e., compund LFPs,
        #CSDs etc.
        populations_path = os.path.join(psc.output_proc_prefix, ps_id,
                                        'populations'),
        
        #location of spike output from the network model
        spike_output_path = os.path.join(psc.output_proc_prefix, ps_id),
        
        #destination of figure file output generated during model execution
        figures_path = os.path.join(psc.output_proc_prefix, ps_id, 'figures')
    ))
    
    
    #parameters for class CachedTopoNetwork instance
    PS.update(dict(network_params = dict(
            simtime = PSET.t_sim,
            dt = PSET.dt,
            spike_output_path = PS.spike_output_path,
            label = PSET.spike_detector_label,
            ext = 'gdf',
            GIDs = PARAMS.GIDs,
            X = PARAMS.X,
            label_positions = PSET.position_filename_trunk,
        )))
    
    
    #population (and cell type) specific parameters
    PS.update(dict(
        #list of presynaptic neuron populations (network populations)
        X = PARAMS.X,
        
        #list of postsynaptic neuron populations
        Y = PARAMS.Y,
        
        #list of cell types y in each postsynaptic population Y 
        y = PARAMS.y,
        
        #population-specific LFPy.Cell parameters
        cellParams = PARAMS.yCellParams,
        
        #assuming excitatory cells are pyramidal
        rand_rot_axis = PARAMS.rand_rot_axis,
        
        #population sizes
        N_X = PARAMS.N_X,
        
        #kwargs passed to LFPy.Cell.simulate()
        simulationParams = dict(),
        
        #set up parameters corresponding to cylindrical model populations
        populationParams = PARAMS.populationParams,
        
        #set the boundaries between the "upper" and "lower" layer
        layerBoundaries = PARAMS.layerBoundaries,
        
        #set the geometry of the virtual recording device                       
        electrodeParams = dict(
                #contact locations:
                x = np.mgrid[-18:19:4, -18:19:4][1].flatten()*100,
                y = np.mgrid[-18:19:4, -18:19:4][0].flatten()*100,
                z = np.array([PARAMS.layerBoundaries.mean(axis=1)[1] for x in range(100)]), #center of layer 2/3
                #extracellular conductivity:
                sigma = 0.3,
                #contact surface normals, radius, n-point averaging
                N = [[0, 0, 1]]*100,
                r = 5,
                n = 50,
                seedvalue = None,
                #dendrite line sources, soma as sphere source (Linden et al 2014),
                #and account for periodic boundary conditions also in the LFP to
                #the second order. Must use LFPy from github.com/espenhgn and
                #branch som_as_point_periodic (as the periodic boundaries are
                #just a hack)
                method = 'soma_as_point_periodic' if PSET.pbc else 'soma_as_point',
                periodic_order=2,
                periodic_side_length=PSET.extent_length*1000.,                
        ),
        
        #runtime, cell-specific attributes and output that will be stored
        savelist = [
            'somav',
            'somapos',
            'x',
            'y',
            'z',
            'LFP',
        ],
        
        #flag for switching on calculation of CSD
        calculateCSD = False,
        
        #time resolution of saved signals
        #TODO: find nicer way to make sure we use same dt in analysis and
        #LFP simulation output
        dt_output = PSET.dt*5 
    ))
    
     
    #for each population, define layer- and population-specific connectivity
    #parameters
    PS.update(dict(
        #number of connections from each presynaptic population onto each
        #layer per postsynaptic population, preserving overall indegree
        k_yXL = PARAMS.k_yXL,
        
        #set up table of synapse weights onto postsynaptic cell type y from each
        #possible presynaptic population X 
        J_yX = PARAMS.J_yX,
        
        
        #table of synapse time constants onto postsynaptic cell type y from each
        #possible presynaptic population X
        tau_yX = PARAMS.tau_yX,
        
        #unit synapse weights used for LFP proxy
        J_yX_unit = PARAMS.J_yX_unit,
    
        #set up synapse parameters as derived from the network
        synParams = PARAMS.synParams,
        
        #set up delays, here using fixed delays of network unless value is None
        synDelayLoc = {y : [None for X in PS.X] for y in PS.y},
        
        #distribution of delays added on top of per connection delay using either
        #fixed or linear distance-dependent delays
        synDelayScale = {y : [PSET.sigma_delays['E' not in X] for X in PS.X] for y in PS.y},
        
        
        #For topology-like connectivity. Only exponential and gaussian
        #connectivity and circular
        #masks are supported with fixed indegree given by k_yXL,
        #using information on extent and edge wrap (periodic boundaries). 
        #At present, synapse delays are not distance-dependent. For speed,
        #multapses are always allowed. 
        topology_connections = PARAMS.topology_connections,
    ))
    
    
    #mapping between population type and cell type specificity,
    PS.update(dict(
        mapping_Yy = PARAMS.mapping_Yy
    ))
    
    
    ###### MISC ATTRIBUTES NOT USED FOR SIMULATION ####
    PS.update(dict(
        depths = PARAMS.depths,
        PATH_m_y = PARAMS.PATH_m_y,
        m_y = PARAMS.m_y,
        N_y = PARAMS.N_y,
        full_scale_num_neurons = PARAMS.full_scale_num_neurons,
        
    ))
    
    #### Attributes for LFP proxy stuff####
    PS.update(dict(
        #CachedFixedSpikesTopoNetwork params
        activationtimes=[200, 300, 400, 500, 600, 700, 800, 900, 1000],
        filelabel='population_spikes',
        mask=[-400, 0, -400, 0],
        #output files for TopoPopulation and Postprocess classes
        output_file = '{}_population_{}',
        compound_file='{}sum.h5',
        #parameters for LFP kernel extraction
        nlag = int(20 / PS.dt_output),
        
        #output files for calculating LFP proxies
        output_file_proxy = '{}_proxy_{}',
        compound_file_proxy='{}proxy.h5',
        #output files for estimates from firing rates
        output_file_approx = '{}_population_approx_{}',
        compound_file_approx='{}_approx.h5',
        extent = PSET.extent_length*1E3,
        pbc = PSET.pbc,
        
        #synapse strength in units of pA
        J = PARAMS.J
    ))
    
    return PS, PSET
