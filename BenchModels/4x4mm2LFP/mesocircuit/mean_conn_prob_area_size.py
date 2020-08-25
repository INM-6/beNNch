#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Calculations based on Sheng, T. K. The Distance between Two Random
Points in Plane Regions, Advances in Applied Probability,
Vol. 17, No. 4 (Dec., 1985), pp. 748-773

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
import scipy
from scipy import integrate
import json
import os


def get_sigma_C_0_PD14():
    '''
    Compute sigma (mm) and C_0 for the Potjans&Diesmann values.
    PD14 Eqs. 7&8
    '''
    ### compute sigma and C_0 from C_p, C_a and r_p ###

    # mean values from Potjans&Diesmann2014, table S1
    C_p = 0.139
    C_a = 0.079
    # sampling radius of physiological map in mm
    r_p = 0.1

    sigma = r_p / np.sqrt(-2.*np.log(1. - (np.pi*r_p**2*C_p) / (C_a) ))

    C_0 = C_a / (2. * np.pi * sigma**2)

    return sigma, C_0


def get_C_mean_PD14(R, sigma, C_0):
    '''
    Compute mean connection probability as done by Potjans&Diesmann.
    PD14 Eq. 9
    '''
    c_mean = 2./R**2 * C_0 * sigma**2 * (1. - np.exp(-R**2 / (2.*sigma**2)))
    return c_mean


# assume circular layers, theorem 2.4
def get_C_mean_circle(R, sigma, C_0):
    '''
    Compute the mean connection probability for a circular layer,
    taking into account the distances between all neuron pairs
    in that domain.

    Parameters
    ----------
        R : float
            radius of the circle in mm
        sigma : float
            lateral spread in mm
        C_0 : float
            connection probability at 0 distance
    Returns
    -------
        C_mean : float
            mean connection probability
    '''

    def integrand(r, R, sigma):
        return r * np.exp(-r**2/(2*sigma**2)) * \
            ( 4*scipy.arctan(np.sqrt( (2*R-r)/(2*R+r) )) - \
              scipy.sin(4*scipy.arctan(np.sqrt( (2*R-r)/(2*R+r) ))) )

    limits = [0., 2.*R]

    C_mean = 2*C_0/(np.pi*R**2) * scipy.integrate.quad(integrand,
                                                       limits[0],
                                                       limits[1],
                                                       args=(R, sigma))[0]
    return C_mean


# assume square layers, theorem 2.6
def get_C_mean_square(L, sigma, C_0):
    '''
    Compute the mean connection probability for a square layer,
    taking into account the distances between all neuron pairs
    in that domain.

    Parameters
    ----------
        L : float
            side length of square in mm
        sigma : float
            lateral spread in mm
        C_0 : float
            connection probability at 0 distance
    Returns
    -------
        C_mean : float
            mean connection probability
    '''

    def integrand1(r, L, sigma):
        return np.exp(-r**2/(2.*sigma**2)) * r * \
            (np.pi * L**2 - 4.*L*r+ r**2)

    limits1 = [0., L]

    I1 = scipy.integrate.quad(integrand1,
                              limits1[0],
                              limits1[1],
                              args=(L, sigma))[0]

    def integrand2(r, L, sigma):
        return np.exp(-r**2/(2.*sigma**2)) * r * \
            (2.*L**2 * (scipy.arcsin(L/r) - scipy.arccos(L/r) - 1.) + \
             4.*L*r * np.sqrt(1.-(L/r)**2) - r**2)

    limits2 = [L, np.Inf]

    I2 = scipy.integrate.quad(integrand2,
                              limits2[0],
                              limits2[1],
                              args=(L, sigma))[0]

    # normalization was missing in the paper (eq. 2.6.5 is wrong),
    # added division by R**4
    C_mean = 2*C_0 / L**4 * (I1 + I2)

    return C_mean


def get_C_mean_PBC(L, sigma, C_0):
    '''
    TODO implement checks, generate if not exist
    '''
    try:
        script_dir = os.path.dirname(__file__)
        file_name = 'c_mean_pbc.json'
        abs_file_path = os.path.join(script_dir, file_name)
        with open(abs_file_path, 'r') as f:
            data = f.read()
        d = json.loads(data)
        key = 'sigma{:.2f}_cpeak{:.2f}'.format(sigma, C_0)
        C_mean = d[key][str(L)]

    except:
        raise Exception('Fix get_C_mean_PBC.')
    return C_mean


def compute_C_mean_PBC(L_list, sigma, C_0):
    '''
    '''
    def c_gauss(x0, x1, x2, x3): # y2, x2, y1, x1
        xx1 = x3
        xx2 = x1
        yy1 = x2
        yy2 = x0
        r = np.sqrt((xx2-xx1)**2 + (yy2-yy1)**2)
        return C_0 * np.exp(-r**2 / (2*sigma**2))

    def c_mean(L):
        def lim0(x1, x2, x3):
            return [x2-L/2., x2+L/2.]
        def lim1(x2, x3):
            return [x3-L/2., x3+L/2.]

        I = 1./L**4 * integrate.nquad(c_gauss, [lim0, lim1, [-L/2., L/2.], [-L/2., L/2.]])[0]
        return I


    key = 'sigma{:.2f}_cpeak{:.2f}'.format(sigma, C_0)
    d = {}
    d[key] = {}
    for L in L_list:
        d[key][L] = c_mean(L)

    with open('c_mean_pbc.json', 'w') as f:
        json.dump(d, f)

    return


def compute_C_mean_area(area, sigma, C_0):
    '''
    Just call C_mean_circle() and C_mean_square() and print results.

    Parameters
    ----------
        area : float
            area size in mm^2
        sigma : float
            lateral spread in mm
        C_0 : float
            connection probability at 0 distance
    Returns
    -------
    '''

    # radius of circle
    R_circle = np.sqrt(area / np.pi)
    # side length of square
    R_square = np.sqrt(area)

    # compute mean connection probabilities
    C_circle = get_C_mean_circle(R_circle, sigma, C_0)
    C_square = get_C_mean_square(R_square, sigma, C_0)

    print 'area (mm^2) = ', area
    print 'R_circle (mm) = ', R_circle
    print 'R_square (mm) = ', R_square
    print 'C_circle = ', C_circle
    print 'C_square = ', C_square
    print ''

    return


def get_mean_delay_circle(area, d0, v, sigma):
    '''
    Compute the mean delay for a given circular area (= disk).
    Delays have a linear distance dependence.
    Gaussian connection probability.

    Paramaters
    ----------
    area : float
        area size in mm^2
    d0 : float
        delay offset in ms
    v : float
        conduction velocity in mm/ms
    sigma : float
        standard deviation of Gaussian profile in mm

    Returns
    -------
    dmean : float
        mean delay in ms
    '''
    # radius of circle
    R = np.sqrt(area / np.pi)

    limits = [0., 2.*R]

    def integrand_delay(r, d0, v, sigma, R):
        atan = 4. * scipy.arctan(np.sqrt((2.*R - r) / (2.*R + r)))
        return (d0 + r/v) * \
            np.exp(-r**2/(2.*sigma**2)) * \
            r * (atan - np.sin(atan))

    def integrand_Cnorm(r, sigma, R):
        atan = 4. * scipy.arctan(np.sqrt((2.*R - r) / (2.*R + r)))
        return np.exp(-r**2/(2.*sigma**2)) * \
            r * (atan - np.sin(atan))

    I_delay = scipy.integrate.quad(integrand_delay,
                                   limits[0],
                                   limits[1],
                                   args=(d0, v, sigma, R))[0]

    # normalization
    I_Cnorm = scipy.integrate.quad(integrand_Cnorm,
                                   limits[0],
                                   limits[1],
                                   args=(sigma, R))[0]

    dmean = I_delay / I_Cnorm

    return dmean
