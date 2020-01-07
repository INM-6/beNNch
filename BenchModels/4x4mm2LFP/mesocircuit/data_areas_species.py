#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Collection of data from the Potjans&Diesmann2014 model.

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
import copy
import numpy as np

# data valid for 1 mm^2 model, assuming to allow multapses
data = {}

datasets = ['potjans'] # more to be added
data['datasets'] = datasets

keys = ['neuron_numbers',
        'neuron_numbers_th',
        'connection_probabilities',
        'connection_probabilities_th',
        'K_bg',
        'weight_mod_facts',
        'bg_rate',
        'PSP_e',
        'PSP_ext',
        'tau_syn']
for k in keys:
    data[k] = {}

# Potjans&Diesmann2014
data['neuron_numbers']['potjans'] = [[20683, 5834], # L23 E, I
                                     [21915, 5479], # L4 E, I
                                     [4850, 1065],  # L5 E, I
                                     [14395, 2948]] # L6 E, I

# connection probabilities from paper
# columns correspond to source populations; rows to target populations
# source      23E    23I    4E      4I      5E      5I      6E      6I
data['connection_probabilities']['potjans'] = \
    [[0.1009,  0.1689,  0.0437,  0.0818,  0.0323,  0.,      0.0076,  0.    ],
     [0.1346,  0.1371,  0.0316,  0.0515,  0.0755,  0.,      0.0042,  0.    ],
     [0.0077,  0.0059,  0.0497,  0.135,   0.0067,  0.0003,  0.0453,  0.    ],
     [0.0691,  0.0029,  0.0794,  0.1597,  0.0033,  0.,      0.1057,  0.    ],
     [0.1004,  0.0622,  0.0505,  0.0057,  0.0831,  0.3726,  0.0204,  0.    ],
     [0.0548,  0.0269,  0.0257,  0.0022,  0.06,    0.3158,  0.0086,  0.    ],
     [0.0156,  0.0066,  0.0211,  0.0166,  0.0572,  0.0197,  0.0396,  0.2252],
     [0.0364,  0.001,   0.0034,  0.0005,  0.0277,  0.008,   0.0658,  0.1443]]

# weight modification factors from paper
data['weight_mod_facts']['potjans'] = \
    [[1.0, -4.0, 2.0, -4.0, 1.0, -4.0, 1.0, -4.0],
     [1.0, -4.0, 1.0, -4.0, 1.0, -4.0, 1.0, -4.0],
     [1.0, -4.0, 1.0, -4.0, 1.0, -4.0, 1.0, -4.0],
     [1.0, -4.0, 1.0, -4.0, 1.0, -4.0, 1.0, -4.0],
     [1.0, -4.0, 1.0, -4.0, 1.0, -4.0, 1.0, -4.0],
     [1.0, -4.0, 1.0, -4.0, 1.0, -4.0, 1.0, -4.0],
     [1.0, -4.0, 1.0, -4.0, 1.0, -4.0, 1.0, -4.0],
     [1.0, -4.0, 1.0, -4.0, 1.0, -4.0, 1.0, -4.0]]

# background rate from paper
data['bg_rate']['potjans'] = 8.

# postsynaptic potential
data['PSP_e']['potjans'] = 0.15

# external postsynaptic potential
data['PSP_ext']['potjans'] = 0.15

# synaptic time constant (same for exc and inh)
data['tau_syn']['potjans'] = 0.5

# number of thalamic neurons
data['neuron_numbers_th']['potjans'] = [[902]]

# connection probabilities thalamus
# (from thalamic neurons to cortical populations)
data['connection_probabilities_th']['potjans'] = \
    [[0.],     # 23E
     [0.],     # 23I
     [0.0983], # 4E
     [0.0619], # 4I
     [0.],     # 5E
     [0.],     # 5I
     [0.0512], # 6E
     [0.0196]] # 6I

# indegrees for background input
data['K_bg']['potjans'] = [[1600,   # 2/3e
                            1500],  # 2/3i
                           [2100,   # 4e
                            1900],  # 4i
                           [2000,   # 5e
                            1900],  # 5i
                           [2900,   # 6e
                            2100]]  # 6i
