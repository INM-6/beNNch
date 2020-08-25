#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Computation of derived network parameters.

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

import numpy as np
import data_areas_species
import mean_conn_prob_area_size as cpas


def get_all_derived_params(paramset):
    '''
    Main function to get all derived parameters.
    Used in parameterspace_control.py

    Arguments
    ---------
    paramset : dict
        parameterset dictionary

    Returns
    -------
    params : dict
        dict with additional parameters
    '''
    params = {}

    # neuron numbers for full network
    params['full_num_neurons'] = get_neuron_numbers(paramset,
                                                    key='neuron_numbers')
    params['full_num_neurons_th'] = get_neuron_numbers(paramset,
                                                       key='neuron_numbers_th')

    # connection_probabilities for full network
    full_conn_probs = get_connection_probs(paramset,
                                           key='connection_probabilities')
    full_conn_probs_th = get_connection_probs(paramset,
                                              key='connection_probabilities_th')


    # compute synapse numbers from connection probabilities and neuron numbers
    params['full_num_synapses'] = \
        conn_probs_to_syn_num(n_pre=params['full_num_neurons'],
                              n_post=params['full_num_neurons'],
                              conn_prob=full_conn_probs,
                              multapses=True)


    params['full_num_synapses_th'] = \
        conn_probs_to_syn_num(n_pre=params['full_num_neurons_th'],
                              n_post=params['full_num_neurons'],
                              conn_prob=full_conn_probs_th,
                              multapses=True)

    # weight modification factors
    params['weight_mod_facts']  = get_weight_mod_facts(paramset=paramset)

    # adjust external indegrees to preserve the mean input of the
    # Potjans&Diesmann microcircuit
    params['K_bg'] = get_K_bg(paramset, params['full_num_neurons'],
                              params['full_num_synapses'],
                              params['weight_mod_facts'])

    # further modify parameters without compensation by changing K_bg,
    # checks if any _no_comp parameters are in paramset,
    # useful for parameter scans
    params = mod_params_without_compensation(paramset, params)

    # set decay_params if decay_param_ex and/or decay_param_in are given
    params['decay_params'] = get_decay_params(paramset)

    return params


def conn_probs_to_syn_num(n_pre, n_post, conn_prob, multapses=True):
    ''' takes a matrix of pairwise connection probabilities and
        converts it to a matrix of synapse numbers
        for 4 layers and 2 populations per layer

        Parameters
        ----------
        n_pre : source neuron numbers as list of floats (e.g., shape (4.2))
        n_post : target neuron numbers as list of floats (e.g., shape (4.2))
        conn_prob : matrix of pairwise connection probabilities
            as list of shape (8,8), containing floats,
            columns=sources, rows=targets
        multapses : whether to allow multapses, True by default
    
        Returns
        -------
        matrix : matrix of synapse numbers between populations
            as list of shape (8,8), containing floats,
            columns=sources, rows=targets
    '''
    matrix = np.empty_like(conn_prob)
    for src in range(matrix.shape[1]): # sources in columns
        [s0], [s1] = np.unravel_index([src], np.shape(n_pre))
        num_src = n_pre[s0][s1]
        for tgt in range(matrix.shape[0]): # targets in rows
            [t0], [t1] = np.unravel_index([tgt], np.shape(n_post))
            num_tgt = n_post[t0][t1]
            if multapses:
                m = np.log(1.-conn_prob[tgt][src])/np.log((num_src*num_tgt-1.)/(num_src*num_tgt))
            else:
                m = conn_prob[tgt][src] * num_src*num_tgt
            matrix[tgt][src] = round(m, 6)
    return matrix.tolist()


def get_neuron_numbers(paramset, key='neuron_numbers'):
    '''
    Adapts neuron numbers to square sheet of side length "extent_length".

    Arguments
    ---------
    paramset : dict
        parameterset dictionary
    key : str
        "neuron_numbers" or "neuron_numbers_th" (for thalamic neurons)

    Returns
    -------
    full_num_neurons : list of list(s)
        neuron numbers for full network
    '''
    mm2_num_neurons = list(np.array(data_areas_species.data[key][paramset['dataset_exp']]))
    full_num_neurons = [ [int(i * paramset['extent_length']**2) \
                          for i in mm2_num_neurons[j] ] \
                         for j in range(len(mm2_num_neurons)) ]

    return full_num_neurons


def get_connection_probs(paramset, key='connection_probabilities'):
    '''
    Computes connection probabilties for the full network on a square area.
    Applies multiplicative modifications to cortical-cortical connections.
    Uses connection probabilities from original microcircuit if
    random connectivity and 1mm^2.

    Arguments
    ---------
    paramset : dict
        parameterset dictionary
    key : str
        "connection_probabilties" or "connection_probabilities_th"
        (for thalamic neurons)

    Returns
    -------
    full_conn_probs : list of list(s)
        connection probabilities for full network
    '''
    # prevent data_areas_species.data[key][paramset['dataset_exp']]
    # from being overwritten
    mm2_C_dataset_exp = list(np.array(data_areas_species.data[key][paramset['dataset_exp']]))

    # use connection probabilities from original microcircuit if
    # random connectivity and 1mm^2
    if paramset['connect_method']=='Connect' and paramset['extent_length']==1:
        full_conn_probs = mm2_C_dataset_exp
    else:
        _sigma, _C_0 = cpas.get_sigma_C_0_PD14() # from PD14 Eqs. 7&8

        # next, compute the mean connection probability for this area
        if paramset['pbc']:
            C_mean = cpas.get_C_mean_PBC(paramset['extent_length'], _sigma, _C_0)
        else:
            C_mean = cpas.get_C_mean_square(paramset['extent_length'], _sigma, _C_0)

        # mean connection probability of the original map, 1mm^2
        # mm2_C_mean_dataset_exp = np.mean(mm2_C_dataset_exp)
        mm2_C_mean  = cpas.get_C_mean_PD14(1./np.sqrt(np.pi), _sigma, _C_0) # from PD14 Eq. 9

        # scale connection probabilities from the original map to this area
        full_conn_probs = [ [ i * C_mean/mm2_C_mean for i in mm2_C_dataset_exp[j] ] \
                            for j in range(len(mm2_C_dataset_exp)) ]

    # linear scaling of connection probabilities
    if key=='connection_probabilities':
        for i, j, corr in paramset['conn_prob_modifications']:
            # modifications assumed to be in the form [post, pre, correction]
            full_conn_probs[i][j] = full_conn_probs[i][j] * corr
    elif key=='connection_probabilities_th':
        conn_prob_factor = paramset['scale_th_conn_prob']
        full_conn_probs = [[c * conn_prob_factor for c in i] for i in full_conn_probs]

    return full_conn_probs


def get_weight_mod_facts(paramset):
    '''
    Creates array of weight modification factors relative to PSP_e
    for all pre- and postsynaptic cortical popoulations.
    IPSP amplitude relative to EPSP amplitude is usually g.
    g=-4 in the microcircuit.
    If dataset_exp is potjans or mam_macaque_V1, the exception to double
    the weight L4e->L2/3e is included as in the original microcircuit.
    Finally, additional weight_modifications are applied.

    Parameters
    ----------
    paramset : dict
        parameterset dictionary

    Returns
    -------
    data : np.ndarray
        2D array with weight modification factors

    '''
    data = np.ones((8, 8), dtype=float)
    data[:, 1::2] *= paramset['g']
    if paramset['dataset_exp'] == 'potjans' or paramset['dataset_exp'] == 'mam_macaque_V1':
        data[0, 2] *= 2
    else:
        raise ValueError('dataset %s not recognized' % paramset['dataset_exp'])

    # TODO: NOT TESTED, YET, ALSO LOOK AT K_BG CHANGES
    # weight change to compensate for change in tau_syn
    delta_g_ex, delta_g_in = get_delta_g_changed_tau_syn(paramset)
    data[:,::2] *= delta_g_ex
    data[:,1::2] *= delta_g_in

    # old version
    #data[:,1::2] *= paramset['tau_syn_ex']/paramset['tau_syn_in']

    for (post, pre, factor) in paramset['weight_modifications']:
        data[post, pre] *= factor

    return data.tolist()


def get_K_bg(paramset, full_num_neurons, full_num_synapses, weight_mod_facts):
    '''
    Computes the numbers of external synapses for preserving the firing rates
    "potjans_mean_rates".
    Applies multiplicative modifications to K_bg.

    Arguments
    ---------
    paramset : dict
        parameterset dictionary
    full_num_neurons : list of lists
        neuron numbers for full network
    full_num_synapses : list of lists
        synapse numbers for full network

    Returns
    -------
    K_bg : list
        numbers of external synapses per population
    '''

    # number of external synapses from dataset for 1 mm^2
    K_bg = data_areas_species.data['K_bg'][paramset['dataset_exp']]

    # number of neurons and synapses from dataset
    mm2_num_neurons = list(np.array(data_areas_species.data['neuron_numbers'][paramset['dataset_exp']]))
    mm2_C_dataset_exp = list(np.array(data_areas_species.data['connection_probabilities'][paramset['dataset_exp']]))
    mm2_num_synapses = conn_probs_to_syn_num(mm2_num_neurons, mm2_num_neurons,
                                             mm2_C_dataset_exp,
                                             multapses=True)

    # weight modification factors from dataset
    mm2_weight_mod_facts = list(np.array(data_areas_species.data['weight_mod_facts'][paramset['dataset_exp']]))

    # background rate from dataset
    mm2_bg_rate = data_areas_species.data['bg_rate'][paramset['dataset_exp']]

    # postsynaptic potential, equal to PSP_ext in Potjans&Diesmann
    mm2_PSP_e = data_areas_species.data['PSP_e'][paramset['dataset_exp']]
    mm2_PSP_ext = data_areas_species.data['PSP_ext'][paramset['dataset_exp']]

    #print np.sum(mm2_num_neurons), np.sum(full_num_neurons)
    #print np.sum(mm2_num_synapses), np.sum(full_num_synapses)

    # change applied to weight_mod_facts based on changed synaptic time constants
    delta_g = get_delta_g_changed_tau_syn(paramset) # [delta_g_ex, delta_g_in]

    # scaling of the number of external synapses dependent on the new indegrees
    # in order to preserve the firing rates
    dim1,dim2 = np.shape(full_num_synapses)
    new_K_bg = np.zeros_like(K_bg).tolist()
    for i in range(dim1): # to

        sum_j = 0.

        for j in range(dim2): # from

            # indegree * weight modification * PSP_e
            full_input = full_num_synapses[i][j] / full_num_neurons[i // 2][i % 2] \
                         * weight_mod_facts[i][j] \
                         / delta_g[j % 2] \
                         * paramset['PSP_e']

            mm2_input = mm2_num_synapses[i][j] / mm2_num_neurons[i // 2][i % 2] \
                        * mm2_weight_mod_facts[i][j] \
                        * mm2_PSP_e

            diff_input = mm2_input - full_input

            # sum differences in recurrent input
            sum_j += diff_input * paramset['potjans_mean_rates'][j // 2][j % 2]


        # assume weight_mod_fact is 1 for external input
        new_K_bg[i // 2][i % 2] = 1./(paramset['PSP_ext'] * paramset['bg_rate']) \
                                  * (K_bg[i // 2][i % 2] * mm2_PSP_ext * mm2_bg_rate \
                                     + sum_j)

    # the K_bg_modifications must be applied after other corrections
    # based on the Potjans rates.
    # apply modifications to noise indegree
    for pop,mod in paramset['K_bg_modifications']:
        i,j = np.unravel_index(pop, (4,2))
        new_K_bg[i][j] *= mod

    return new_K_bg


def get_delta_g_changed_tau_syn(paramset):
    '''
    Changes to be applied to weight_mod_facts
    based on changed synaptic time constants.
    '''
    def term(tm, ts):
        return (tm/ts)**(tm/(ts-tm)) - (tm/ts)**(ts/(ts-tm))
    def delta_g(tm, ts0, ts1):
        return (ts0-tm)/(ts1-tm) * term(tm, ts1)/term(tm, ts0)

    tau_syn0 = data_areas_species.data['tau_syn'][paramset['dataset_exp']]
    delta_g_ex = delta_g(paramset['tau_m'], tau_syn0, paramset['tau_syn_ex'])
    delta_g_in = delta_g(paramset['tau_m'], tau_syn0, paramset['tau_syn_in'])

    return [delta_g_ex, delta_g_in]


def mod_params_without_compensation(paramset, params):
    '''
    Checks if parameters bg_rate_comp_nocomp and g_comp_nocomp exist in
    parameterspace. They are lists containing the original (already compensated
    for by adapting K_bg) and an updated value. The original parameter
    is then overwritten by the updated one.
    '''
    if 'bg_rate_comp_nocomp' in paramset:
        if paramset['bg_rate_comp_nocomp'][0] != paramset['bg_rate']:
            raise Exception('First entry of bg_rate_comp_nocomp '
                            +'must be bg_rate for which we compensate.')

        # overwrite bg_rate with the new uncompensated rate
        params['bg_rate'] = paramset['bg_rate_comp_nocomp'][1]

    if 'g_comp_nocomp' in paramset:
        if paramset['g_comp_nocomp'][0] != paramset['g']:
            raise Exception('First entry of g_comp_nocomp '
                            +'must be g for which we compensate.')

        # overwrite weight_mod_facts and g
        wmf = np.array(params['weight_mod_facts'])
        wmf[:, 1::2] *= paramset['g_comp_nocomp'][1] / paramset['g']
        params['weight_mod_facts'] = wmf.tolist()
        params['g'] = paramset['g_comp_nocomp'][1]

    return params


def get_decay_params(paramset):
    '''
    If decay_param_ex/in in paramset, matrix decay_params is overwritten
    with specified values.
    '''
    decay_params = np.array(paramset['decay_params'])
    if 'decay_param_ex' in paramset and paramset['decay_param_ex']!=None:
        decay_params[:, ::2] = np.tile([paramset['decay_param_ex']],
                                       (len(decay_params), len(decay_params)/2))
    if 'decay_param_in' in paramset and paramset['decay_param_in']!=None:
        decay_params[:, 1::2] = np.tile([paramset['decay_param_in']],
                                       (len(decay_params), len(decay_params)/2))

    return decay_params.tolist()
