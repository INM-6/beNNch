#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Parameter definition for mesocircuit model.

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
import numpy as np
from NeuroTools import parameters as ps

# independent parameters
base_parameters = ps.ParameterSet(dict(

    ###################################################
    ###             Simulation parameters           ###
    ###################################################
    #master seed for random number generators
    #actual seeds will be master_seed ... master_seed + 2*n_vp
    # ==>> different master seeds must be spaced by at least 2*n_vp + 1
    #see Gewaltig et al (2012) for details
    master_seed = 123456, # changes rng_seeds and grng_seed
    #number of virtual processes, must be even dividable by # MPI SIZE
    #n_vp = ps.ParameterRange([96]), #even dividable with 16 and 24 cores
    #n_vp = lnodes*ppn, #number of virtual processes

    # full simulation time in NEST is transient + t_sim, but from the
    # time assigned to the transient, no spikes are recorded
    # simulation time (ms)
    t_sim  = 500.,#2000.,
    # transient (ms)
    transient = 50.,#500.,

    dt = 0.1, # simulation step (ms). Default is 0.1 ms.


    ####################################################
    ####             Network parameters              ###
    ####################################################

    # reduce the neuron density, e.g., for test simulations
    density = 0.005,#1.,

    # layers and populations
    layers = ['L23', 'L4', 'L5', 'L6'],
    populations = ['e', 'i'],

    # use neuron numbers and connectivity data from a dataset in
    # data_areas_species.py (currently, choose between 'potjans',
    # 'mam_macaque_V1' and 'chadderdon')
    dataset_exp = 'potjans',

    # neuron model.
    # For PSP-to-PSC conversion to be correct, synapses should be
    # current-based with an exponential time course
    neuron_model = '/iaf_psc_exp',

    Vm0_mean = -58.0, # mean of initial membrane potential (mV)
    Vm0_std = 10.0, # std of initial membrane potential (mV)

    # neuron model parameters
    # The model parameter dictionary is now set up in parameterspace_control.py.
    # I want to be able to vary certain parameters in the parameter space,
    # however, these are also at times model neuron specific like tau_syn_ex.
    # Thus, each model neuron we want to use must have its own set of
    # parameters below.
    tau_m = 10.,        # membrane time constant (ms)
    tau_syn_ex = 0.5,   #excitatory synaptic time constant (ms)
    tau_syn_in = 0.5,   #inhibitory synaptic time constant (ms)
    t_ref = 2.,         # absolute refractory period (ms)
    E_L = -65.,         # resting membrane potential (mV)
    V_th = -50.,        # spike threshold (mV)
    C_m = 250.,         # membrane capacitance (pF)
    V_reset = -65.,     # reset potential (mV)

    # additional model parameters, here we just use this as a container (perhaps
    # filled with defaults).
    # model_params = {},


    # type of synapse
    synapse_model = '/static_synapse',

    # weights

    # mean EPSP amplitude (mV) of connections before multiplication with
    # weight_mod_facts (e.g., weight_mod_facts != 1 will modify the weight)
    PSP_e = 0.15,

    # global factor for inhibitory scaling,
    # changes are compensated by adapting K_bg.
    # For compensating according to g and then using another inhibitory scaling,
    # use g_no_comp in addition.
    g = -4.,

    # Modify connection probabilities of intracortical connections while
    # attempting to retain working point similar to Bos et al., 2016, by
    # modifying noise indegrees. These conn_prob_modifications will be
    # multiplied with the corresponding connection probability specified in the
    # file data_areas_species.py.
    # An entry "[3, 2, 0.9]" (i.e., [post, pre, correction]) would mean that
    # connection probabilities
    # of connections to L4I from L4E is reduced by 10 %.
    # The mean input is preserved by adapting the noise indegree to L4I.
    conn_prob_modifications = [],

    #modification of per population weights. To for example half the weight
    #of L2/3E -> L4I, the entry would be [[3, 0, 0.5]], in short
    #[[i_post, i_pre, factor]]
    weight_modifications = [],

    #modification of per population noise in-degree, as the default noise
    #indegrees are defined in data_areas_species. These float factors will be
    #multiplied with the default noise indegree and converted to int.
    # An entry "[3, 0.9]" (i.e., [population, correction]) means that
    # the noise indegree of population 3, that is L4I, is reduced by 10 %.
    # K_bg_modifications are applied after automatically
    # preserving the mean input accounting for changes in PSP_e,
    # connection probabilities or weights.
    K_bg_modifications = [],


    # (normal) weights as in the mircorcircuit or (lognormal) ones
    weight_parameters = [ 'normal', 'normal' ],

    # standard deviation of PSC amplitudes relative to mean PSC amplitudes
    # for normal weights
    PSC_rel_sd = [ 0.1, 0.1 ], # exc. and inh.

    # standard deviation of PSC amplitudes relative to mean PSC amplitudes
    # for lognormal weights
    PSC_lognorm_rel_sd = [ 0.051, 0.013 ],


    # dendritic delays for excitatory and inhibitory transmission

    # delay parameter, random (normal) or distance-dependent (linear)
    delay_parameter = 'linear',

    # set normal delays to mean value of linear delays on a disk
    # of area 1mm^2, assuming decay_param=0.3 for all connections.
    # if True and delay_parameter = 'normal',
    # variable 'delays' is overwritten and used.
    # if False and delay_parameter = 'normal',
    # the here entered 'delays' are used.
    use_mean_delays = False,

    # 1. parameters for delay_parameter == 'normal'
    # mean delay (ms)
    delays =  [1.5, 0.75],
    # standard deviation relative to mean delays
    delay_rel_sd = 0.5,

    # 2. parameters for delay_parameter == 'linear'
    # minimum for linear delay for excitatory and inhibitory connections
    # (we used to use 0.3 for both)
    min_delays = [ 0.5, 0.5 ],      # exc. and inh.

    # propagation speed needed for calculating distance-dependent delays
    # ( 1 m/s = 1 mm/ms )
    # computed for preserving the mean delays of the microcircuit
    prop_speeds = [ 0.3, 0.3],     # exc. and inh.

    # sigma for getting variability in delays from normal distribution
    sigma_delays = [ 0.1, 0.1 ],    # exc. and inh.

    # mean rates in the full-scale model, necessary for scaling
    # precise values differ somewhat between network realizations
    potjans_mean_rates = [ [0.7089913,  2.74821237], # microcircuit (theo)
                           [4.56197259, 5.78825561],
                           [7.27856032, 8.46814331],
                           [1.06378518, 7.65750536] ],

    ###################################################
    ###           Stimulus parameters               ###
    ###################################################


    # DC amplitude at each external input synapse (pA)
    # This is relevant for reproducing Potjans & Diesmann (2012) Fig. 7.
    # Usually, it is set to 0.
    dc_amplitudes = [ [0., 0.],
                      [0., 0.],
                      [0., 0.],
                      [0., 0.] ],

    # rate of background Poisson input at each external input synapse
    # (spikes/s),
    # changes are compensated by adapting K_bg.
    # For compensating according to bg_rate and then using another background
    # rate, use in addition bg_rate_comp_nocomp = [bg_rate, new_bg_rate].
    bg_rate = 8.,


    # mean EPSP amplitude (mV) for external input
    # (thalamus + external Poisson generators)
    PSP_ext = 0.15,

    # scaling factor for thalamic connections
    scale_th_conn_prob = 1.,

    # connection probabilities for thalamic input
    conn_profiles_th = [ ['gaussian',    # 2/3e
                          'gaussian'],   # 2/3i
                         ['gaussian',    # 4e
                          'gaussian'],   # 4i
                         ['gaussian',    # 5e
                          'gaussian'],   # 5i
                         ['gaussian',    # 6e
                          'gaussian'] ], # 6i

    # decay parameters for thalamic input
    decay_params_th = [ [0.3,    # 2/3e
                         0.3],   # 2/3i
                        [0.3,    # 4e
                         0.3],   # 4i
                        [0.3,    # 5e
                         0.3],   # 5i
                        [0.3,    # 6e
                         0.3] ], # 6i

    # mask radii for thalamic input
    mask_radii_th = [ [2.,    # 2/3e
                       2.],   # 2/3i
                      [2.,    # 4e
                       2.],   # 4i
                      [2.,    # 5e
                       2.],   # 5i
                      [2.,    # 6e
                       2.] ], # 6i

    # factor for scaling the thalamic weights by multiplication with PSC_ext
    th_weight_scale = 1.,

    # only thalamic neurons within a circle of th_radius are stimulated
    th_radius = 1.,

    # if th_poisson_input, the thalamic popultion is activated at th_start
    # for th_duration and fires with a Poisson rate of th_rate
    # else: all thalamic neurons fire simultaneously at specific th_spike_times
    th_poisson_input = True,
    th_start = 1e5,#100.,     # onset of thalamic input (ms) after transient
    th_duration = 50.,   # duration of thalamic input (ms)
    th_rate = 50.,       # rate of thalamic neurons (spikes/s), Potjans: 15 Hz
    th_spike_times = list(np.arange(100., 1e5, 100.)), # thalamic spike times


    # additional external stimulus with a spike generator to a
    # given subset of neurons of specified populations
    stim_pop = [[0, 0], [1, 0], [0, 0], [0, 0]],
    # center position of spike generator
    stim_pos = [[0., 0.]],
    # mask for selecting a region to be stimulated
    stim_mask = {'/rectangular': {'/lower_left': [-0.5, -0.5],
                                  '/upper_right': [0.5, 0.5]}},
    # factor for scaling the weight of the additional stimulus by
    # multiplication with PSC_ext
    stim_weight_scale = 1.E6,
    # spike times of additional stimulus (ms) after transient
    stim_spike_times = [1.E6],

    ###################################################
    ###       Topology specific parameters          ###
    ###################################################

    # connections routine from NEST: 'Connect' or 'ConnectLayers' (topology)
    connect_method = 'ConnectLayers',

    # completely random positions or alternatively a jittered lattice
    random_positions = True,

    # parameter for neuron positions on jittered lattice (not random)
    # jitter_factor = 0: regular grid points
    # = 1: maximal randomness (max. displacement: +- extent_length/2);
    # random_positions = true and jitter_factor = 1 is comparable to
    # random_positions true,
    # but the neuron numbers are constrained to square numbers in the
    # first case
    jitter_factor = 0., # in [0,1]

    # creation and connection of neurons with topology
    # length of the quadratic topology layer (mm)
    extent_length = 4.,

    # periodic boundary conditions
    pbc = True,

    # choose (uniform), (gaussian) or (exponential) as LOCAL connectivity
    # profile
    # or (random) for position-independent homogeneous connectivity
    #
    # uniform kernel:
    # p(d) = constant
    # gaussian kernel:
    # p(d) = c + p_center*exp(-(d-mean)^2/(2*sigma^2)) with c=mean=0
    # exponential kernel:
    # p(d) = c + a*exp(-d/tau) with c=0
    # /p_center and /a will be set in mesocircuit.sli to preserve the degree
    conn_profiles = np.tile(['gaussian', 'gaussian'], (8,4)).tolist(), # (EX, IN)

    # decay parameters for connection probability dependent on connectivity profile
    # equivalent to tau in exponential profile and sigma in gaussian profile
    decay_params = np.tile([0.3, 0.1], (8,4)).tolist(), # (EX, IN)
    # if decay_param_ex/in are not None, decay_params will be overwritten respectively
    decay_param_ex = 0.3,
    decay_param_in = 0.1,

    # if spatial profile == 'default': conn_profiles and decay_params are
    # specified as above. otherwise, a spatial profile is seleced from
    # data_areas_species.py
    spatial_profile = 'default',

    # set the mask size (cutoff radius) in (mm) (without pbc adjust scale_conn_prob accordingly)
    mask_radii = np.tile([2., 2.], (8,4)).tolist(), # (EX, IN)


    ###################################################
    ###         Recording parameters                ###
    ###################################################

    # whether to record voltage from a fixed fraction of neurons in
    # each population
    overwrite_existing_files = True,

    # whether to record spikes from a fixed fraction of neurons in each
    # population. If false, a fixed number of neurons is recorded in each
    # population.
    # record_fraction_neurons_spikes true with f_rec_spikes 1. records all
    # spikes
    record_fraction_neurons_spikes = True,
    frac_rec_spikes = 1.,
    n_rec_spikes = 2000,

    # whether to record voltage from a fixed fraction of neurons in each
    # population
    # the number of neurons to record voltages from must be a subset of
    # the neurons of which spikes are recorded because we save only the
    # positions of these neurons
    record_fraction_neurons_voltage = True,
    frac_rec_voltage = 0.1,
    n_rec_voltage = 50,
    # start and stop time for recording voltages (ms) after transient
    rec_voltage_start = 700.,
    rec_voltage_stop = 750.,
    dt_voltage = 1., # time interval in ms for recording voltages

    # check the number of desired and established connections
    # between populations
    check_connections_true = False,
    # same for thalamic population
    check_th_connections_true = False,

    # whether to write any recorded cortical spikes to file
    save_cortical_spikes = True,

    # whether to write any recorded membrane potentials to file
    save_voltages = False,

    # whether to record thalamic spikes (only used when n_thal in
    # network_params.sli is nonzero)
    record_thalamic_spikes = True,

    # whether to write any recorded thalamic spikes to file
    save_thalamic_spikes = True,

    # wheter to write any connections to file
    save_connections =  False,

    # number of targets for which connections shall be stored,
    # per rank and per population
    max_num_targets_to_write = 2,
    # maximum number of connections per target
    max_num_conn_to_write = 10000000,

    # name of file to which to write global IDs
    GID_filename = 'population_GIDs.dat',

    # name of file to which to write the neuron positions.
    # in mesocircuit.sli, the rank is appended to get position_filename
    position_filename_trunk = 'neuron_positions_rank',

    # name of file to which to write the connectivity information
    # in mesocircuit.sli, the rank is appended to get connection_filename
    connection_filename_trunk = 'neuron_connections_rank',

    # name of file for writing some relevant simulation times
    sim_info_filename = 'sim_info.dat',

    # stem for spike detector file labels
    spike_detector_label = 'spikes_',

    # stem for voltmeter file labels
    voltmeter_label = 'voltages_',

    # stem for thalamic spike detector file labels
    th_spike_detector_label = 'spikes_8',

    # name of file to which to write derived parameters
    derived_params_filename = 'derived_params.dat',

    # file name for standard output
    std_out = 'output.txt',

    # file name for error output
    error_out = 'errors.txt',

))

ParameterSpaces = {}

# base parameters
ParameterSpaces['a'] = ps.ParameterSpace(base_parameters)




if __name__ == '__main__':
    pass
