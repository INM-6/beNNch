"""
Python version of population model

Population model for NEST 3.0

The population model creates scale * 20 populations, where each population
consist of 5000 neurons. When connecting, we connect each population with
100 other random populations, and use a indegree of 50 connections per neuron.

This model focuses on connection, no simulation is carried out, except for a
short simulation to carry out the exchange needed because of connection.

- Connections use stdp_pl_synapse_hom_hpc if plastic is true, else static_synapse

- Network size: scale * 20 * 5000
- Connectivity: 50 * 100 * num_populations

User arguments
--------------
- All are required
- Script must be run as

  python <script_path>/population_blocks.py scale num_vp plastic

- Arguments:
  - scale (int)    : network size is scale * 20 * 5000 (scale >= 5 recommended)
  - num_vp (int)   : number of virtual processes (MPI ranks x threads per rank)
  - plastic (bool) : if True, use plastic synapses, otherwise static

Benchmark phases
----------------
1. Startup
2. Neuron (node) creation
3. Network connection setup
4. Initial simulation (single time step)
5. Finalization
"""

import math
import nest
import random
import sys
import time


M_INFO = 10

########################### PARAMETER SECTION #################################

if len(sys.argv) != 5:
    raise ValueError("user arguments should be scale, num_vp, plastic and rule")

user_scale = int(sys.argv[1])
user_nvp = int(sys.argv[2])
user_plastic = bool(sys.argv[3])
user_rule = str(sys.argv[4])

print('GIT: ({}) \nuser_scale: {:>3} \nuser_nvp: {:>5} \nuser_plastic: {} \nuser_rule: {:>5}'.format(
    nest.version(), user_scale, user_nvp, user_plastic, user_rule))

if user_rule not in ['all', 'in', 'out', 'tot', 'bern']:
    raise ValueError("User rule not valid.")

# define all relevant parameters: changes should be made here
params = {
    'nvp': user_nvp,      # total number of virtual processes
    'scale': user_scale,  # scaling factor of the network size,
                          # total network size = scale*5000*20 neurons
    'plastic': user_plastic,
    'd_min': 1.5,
    'd_max': 1.5,

    'rule': user_rule,

    'inisimtime': 10., # initial simulation time given in ms: calibration etc
    'dt': 0.1,         # simulation step

    'tau_syn': 0.32582722403722841 # Rise time of synaptic currents in ms
}

logger_params = {'num_nodes': 0}  # Total size of all populations

# Parameter dependencies
brunel_params = {
    'num_neurons': 5000,                  # number of neurons per population
    'num_pop': 20 * params['scale'],      # number of populations (20 * scale)
    'num_pop_connections': 100,           # each population connect to 100 populations
    'model_params': {                     # variables for iaf_psc_alpha
        'E_L':     0.0,                   # Resting membrane potential (mV)
        'C_m':   250.0,                   # Capacity of the membrane (pF)
        'tau_m':  10.0,                   # Membrane time constant (ms)
        't_ref':   0.5,                   # duration of refractory period (ms)
        'V_th':   20.0,                   # Threshold (mV)
        'V_reset': 0.0,                   # Reset Potential (mV)
        'tau_syn_ex': params['tau_syn'],  # time const. postsynaptic currents (ms)
        'tau_minus':  30.,                # time constant for STDP (depression) (ms)
        'V_m': 5.7                        # mean value of membrane potential (mV)
    },
    'delay':  params['d_min'],            # synaptic delay, all connections (ms)
    'stdp_params': {
        'delay':    params['d_min'],
        'alpha':    0.0513,
        'lambda':   0.1,                  # STDP step size
        'mu':       0.4,                  # STDP weight dependence exponent (potentiation)
        'tau_plus': 15.0,                 # time constant for potentiation
    }

}

if params['scale'] < 5:
    # if scale is less than 5, we only connect to num_pop per population
    brunel_params['num_pop_connections'] = 20 * params['scale']
if params['rule'] == 'all':
    brunel_params['num_neurons'] = round(math.sqrt(5000))  # This is wrong
    #brunel_params['num_pop'] = round(math.sqrt(20 * params['scale']))
    #brunel_params['num_pop_connections'] = round(math.sqrt(20 * params['scale']))


############################ FUNCTION SECTION ##################################

def BuildNetwork(logger):
    tic = time.time() # start timer on construction
    
    # set global kernel parameters
    nest.SetKernelStatus({'total_num_virtual_procs': params['nvp'],
                          'resolution': params['dt']})

    master_seed = 101
    n_vp = nest.GetKernelStatus(['total_num_virtual_procs'])[0]

    nest.ll_api.sli_run('<< /grng rngdict/MT19937 :: 101 CreateRNG >> SetKernelStatus')

    rng_seeds = list(range(
        master_seed + 1 + n_vp,
        master_seed + 1 + (2 * n_vp)
        ))
    grng_seed = master_seed + n_vp
    kernel_dict = {
        'grng_seed': grng_seed,
        'rng_seeds': rng_seeds
        }
    nest.SetKernelStatus(kernel_dict)

    nest.SetDefaults('iaf_psc_alpha', brunel_params['model_params'])

    # ------------------------- Build nodes ------------------------------------

    nest.message(M_INFO, 'build_network', 'Creating populations.')

    population_list = [nest.Create('iaf_psc_alpha', brunel_params['num_neurons'])
                       for _ in range(brunel_params['num_pop'])]

    for gc in population_list:
        logger_params['num_nodes'] += len(gc)

    BuildNodeTime = time.time() - tic

    logger.log('{} # build_time_nodes'.format(BuildNodeTime))
    logger.log('{} # virt_mem_after_nodes'.format(memory_thisjob()))

    #-------------------------- Connection -------------------------------------

    tic = time.time()  # start timer for selecting targets

    nest.message(M_INFO, 'build_network', 'Finding target populations.')

    targets = [random.sample(population_list,
                             brunel_params['num_pop_connections'])
               for _ in range(brunel_params['num_pop'])]
    
    FindTargetsTime = time.time() - tic

    logger.log('{} # find_targets_time'.format(FindTargetsTime))

    tic = time.time()  # Start timer for connection time

    conn_degree = 50   # number of connections per neuron
    conn_dict = {'allow_autapses': False, 'allow_multapses': True}

    if params['rule'] == 'in':
        conn_dict.update({'rule': 'fixed_indegree', 'indegree': conn_degree})
    elif params['rule'] == 'out':
        conn_dict.update({'rule': 'fixed_outdegree', 'outdegree': conn_degree})
    elif params['rule'] == 'tot':
        conn_dict.update({'rule': 'fixed_total_number', 'N': conn_degree * brunel_params['num_neurons']})
    elif params['rule'] == 'bern':
        conn_dict.update({'rule': 'pairwise_bernoulli', 'p': conn_degree / brunel_params['num_neurons']})
    elif params['rule'] == 'all':
        conn_dict.update({'rule': 'all_to_all'})

    # Create custom synapse types with appropriate values for our connections
    nest.SetDefaults('static_synapse_hpc', {'delay': brunel_params['delay']})

    if params['plastic']:
        brunel_params['stdp_params'].update({'weight': 1.})
        nest.CopyModel('stdp_pl_synapse_hom_hpc', 'syn_ex_ex', brunel_params['stdp_params'])
    else:
        nest.CopyModel('static_synapse_hpc', 'syn_ex_ex', {'weight': 1.})

    nest.message(M_INFO, 'build_network', 'Connecting populations.')

    if params['d_min'] != params['d_max']:
        delays = {'distribution': 'uniform', 'low': d_min, 'high': d_max}
    else:
        delays = params['d_min']

    for source, target_vec in zip(population_list, targets):
        for target in target_vec:
            nest.Connect(source, target, conn_dict, {'synapse_model': 'syn_ex_ex'})

    # read out time used for building
    BuildEdgeTime = time.time() - tic

    logger.log('{} # build_edge_time'.format(BuildEdgeTime))
    logger.log('{} # virt_mem_after_edges'.format(memory_thisjob()))


def RunSimulation():

    nest.set_verbosity(M_INFO)
    logger = Logger()

    logger.log('{} # virt_mem_0'.format(memory_thisjob()))

    # ----------------------- Network Construction -----------------------------

    BuildNetwork(logger)

    # ---------------- Initial simulation: rig and calibrate -------------------

    tic = time.time()

    nest.Prepare()
    nest.Run(params['inisimtime'])

    InitializationTime = time.time() - tic

    logger.log('{} # init_time'.format(InitializationTime))
    logger.log('{} # virt_mem_after_init'.format(memory_thisjob()))

    # ----------------------- Cleanup and output -------------------------------

    nest.Cleanup()

    logger.log('{} # num_neurons'.format(logger_params['num_nodes']))
    logger.log('{} # num_connections'.format(nest.GetKernelStatus('num_connections')))
    logger.log('{} # min_delay'.format(nest.GetKernelStatus('min_delay')))
    logger.log('{} # max_delay'.format(nest.GetKernelStatus('max_delay')))


def memory_thisjob():
    """Obtain current memory usage"""
    return nest.ll_api.sli_func('memory_thisjob')


class Logger(object):
    """
    logger class used to log memory and timing information from
    network simulations.
    """

    def __init__(self):
        # Nothing to do as we write to stdout
        pass

    def log(self, value):
        print(str(nest.Rank()) + ' ' + value)


nest.ResetKernel()
RunSimulation()
