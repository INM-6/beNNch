#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
General helper functions.

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

# global imports
import numpy as np
import sys
import os
import glob
import copy
import operator
import pickle
import hashlib
import scipy.signal as ss
import scipy.sparse as sp
import h5py
from NeuroTools import parameters as ps


#######################################
### OS COMMANDS                     ###
#######################################

def pwd():
    return os.getcwd()


def mkdir(dir_name):
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)


def ls_files(file_stem):
    return glob.glob(file_stem + '*')


#######################################
### DATA I/O                        ###
#######################################


def dump_sli_params(fname, paramset):
    '''
    Write the paramset dictionary into a sli file.
    
    Arguments
    ---------
    fname : str
        filename
    paramset : dict
        parameter dictionary
    '''
    f = open(fname, 'w')
    for key, value in list(paramset.items()):
        f.write(py2sli(key, value) + '\n')
    f.close()


def get_unique_id(paramset):
    '''
    create a unique hash key for input dictionary
    
    Arguments
    ---------
    paramset : dict
        parameter dictionary
    
    Returns
    -------
    key : str
        hash key
    
    '''
    sorted_params = sort_deep_dict(paramset)
    string = pickle.dumps(sorted_params)
    key = hashlib.md5(string).hexdigest()
    return key


def write_params(fname_txt, paramset):
    '''
    write parameters to nice text file
    
    Arguments
    ---------
    fname_txt : str
        file name
    paramset : dict
        parameter dictionary
    '''
    sorted_params = sort_deep_dict(paramset)
    f = open(fname_txt, 'w')
    for key, value in sorted_params:
        f.write(str(key) + ': ' + str(value) + '\n')
    f.close()


def dump_sparse_to_h5(X, f, data, compression='gzip', compression_opts=2):
    '''
    write sparse matrix entry to hdf5 file under groupname X
    
    Arguments
    ---------
    X : str
        top-level group name
    f : file
        <HDF5 file "filename.h5" (mode r+)>
    data : scipy.sparse.coo.coo_matrix
        <NxM sparse matrix>
    compression : str
        compression strategy, see h5py.File.create_dataset
    compression_opts : int
        compression settings, see h5py.File.create_dataset
    
    '''
    if type(data) == sp.coo_matrix:
        x = data
    else:
        x = data.tocoo()

    group = f.create_group(X)
    dset = group.create_dataset('data_row_col',
                               data=np.c_[x.data, x.row, x.col],
                               compression=compression,
                               compression_opts=compression_opts,
                               maxshape = (None, None))
    dset = group.create_dataset('shape', data=x.shape, maxshape= (None,))
    

def load_h5_to_sparse(X, f):
    '''load sparse matrix stored on COOrdinate format from HDF5.
    
    Arguments
    ---------        
    X : str
        group name, group must contain datasets:
            'data', 'row', 'col' vectors of equal length
            'shape' : shape of array tuple
    f : file
        <HDF5 file "filename.h5" (mode r+)>
    
    Returns
    -------
    data : scipy.sparse.csr.csr_matrix
        
    '''
    data = sp.coo_matrix((f[X]['data_row_col'].value[:, 0],
                          (f[X]['data_row_col'].value[:, 1],
                           f[X]['data_row_col'].value[:, 2])),
                         shape=f[X]['shape'].value) 

    return data.tocsr()


def read_gdf(fname):
    gdf_file = open(fname, 'r')
    gdf = []
    for l in gdf_file:
        data = l.split()
        gdf += [data]
    gdf = np.array(gdf, dtype=object)
    # check data type (float or int) in each column
    # (if mixed, cast to float)
    if gdf.size > 0:
        for col in range(gdf.shape[1]):
            if any ('.' in s for s in gdf[:,col]):
                gdf[:,col] = gdf[:,col].astype(float)
            else:
                gdf[:,col] = gdf[:,col].astype(int)
    return np.array(gdf)
        

def write_gdf(gdf, fname):
    gdf_file = open(fname,'w')
    for line in gdf:
        for i in np.arange(len(line)):
            gdf_file.write(str(line[i]) + '\t')
        gdf_file.write('\n')
    return None


def load_h5_data(path= '',data_type='LFP', y=None, electrode=None, warmup=0., scaling=1.):
    '''
    Types: 'CSD' , 'LFP', 'CSDsum', 'LFPsum'
    '''
    assert y is not None or electrode is not None
    if y is not None:
        f = h5py.File(os.path.join(path, '%s_%ss.h5' %(y,data_type)))
        data = f['data'].value[:,:, warmup:]
        if scaling != 1.:
            np.random.shuffle(data)
            num_cells = int(len(data)*scaling)
            data = data[:num_cells,:, warmup:]

    else:
        f = h5py.File(os.path.join(path, '%ssum.h5' %data_type))
        data = f['data'].value[:, warmup:]

    return data

        

#######################################
### Python to SLI conversion        ###
#######################################


def get_sli_variable(fname, variable):
    '''
    open file fname and search through lines until variable is encountered,
    returning the value of the variable encountered.

    Assumes that sli lines are of the form "/variable value def"

    Arguments
    ---------
    fname : str
        path to sli file
    variable : str
        name of variable (with or without leading "/")

    Returns
    -------
    value : str,
        should be converted to appropriate type, e.g., bool, float etc.
    '''
    index = 1
    #open file
    f = open(fname, 'r')
    for line in f:
        arg = line.split()[0]
        if arg==variable or arg=='/'+variable: # with or without '/' in variable name
            value_line = line.split()
    #close file
    f.close()

    # assume format: (/)variable value def, get out value
    value = ' '.join(value_line[1:-1])

    py_value = _convert_value_sli2py(value)

    return py_value


def py2sli(var_name, py_value):
    '''Converter for parameter definitions from Python to SLI.

    Keyword arguments:
    var_name -- variable name as string
    py_value -- value assigned to the variable in Python
    '''
    sli_out = _convert_value_py2sli(py_value)
    return '/' + var_name + ' ' +  str(sli_out) + ' def'


def _convert_value_py2sli(value):
    '''Converter for values from Python to SLI.
    
    Keyword arguments:
    value -- the value can be numeric, string, boolean, list, or dictionary
    '''
    # numeric
    if type(value) in [float, np.float64, int]:
        sli_value = value

    # string
    elif type(value) == str: #deal with paths starting on root, e.g., /work
        if os.path.isdir(value):
            sli_value = '({0})'.format(value)
        elif value[0] == '/': # model or model parameter
            sli_value = value
        else:
            sli_value = '({0})'.format(value)

    # boolean
    elif type(value) == bool:
        sli_value = str(value).lower()

    # list
    elif type(value) == list:
        if len(value) == 0:
            sli_value = '[]'
        else:
            sli_value = '['
            for val in value:
                sli_value += str(_convert_value_py2sli(val)) + ' '
            sli_value = sli_value[:-1] + ']' # remove last ' '

    # dictionary
    elif type(value) == dict or type(value) == ps.ParameterSet:
        sli_value = '<< '
        for key, val in list(value.items()):
            if key[0] == '/':
                sli_value += key + ' ' + str(_convert_value_py2sli(val)) + ' '
            else:
                sli_value += '/' + key + ' ' + str(_convert_value_py2sli(val)) + ' '
        sli_value = sli_value[:-1] + ' >>' # remove last ' '
    else:
        print(value)
        sys.exit('py2sli does not convert the type of the input argument: ' + str(type(value)))
    return sli_value


def _convert_value_sli2py(value):
    '''
    Converts values from SLI to Python.

    Arguments
    ---------
    value : SLI value as a string
        numeric, string, boolean or list

    Returns
    -------
    py_value : Python value
    '''
    el0 = value[0] # first element

    # numeric
    try:
        int0 = int(el0) # test if numeric
        numeric = True
    except ValueError:
        numeric = False

    if numeric:
        py_value = eval(value)

    # string
    elif el0=='(':
        py_value = value[1:-1] # removing ()

    # boolean
    elif value=='true':
        py_value = True
    elif value=='false':
        py_value = False

    # list
    elif el0=='[':
        py_value = eval(value.replace(' ', ','))

    else:
        print(value)
        sys.exit('sli2py does not convert the type of the input argument')

    return py_value



#######################################
### GENERAL                         ###
#######################################  


def sort_deep_dict(d):
    '''
    sort arbitrarily deep dictionaries into tuples
    
    Arguments
    ---------
    d : dict
    
    Returns:
    x : list of tuples of tuples of tuples ...
    '''
    x = sorted(iter(d.items()), key=operator.itemgetter(0))
    for i, (key, value) in enumerate(x):
        if type(value) == dict or type(value) == ps.ParameterSet:
            y = sorted(iter(value.items()), key=operator.itemgetter(0))
            x[i] = (key, y)
            for j, (k, v) in enumerate(y):
                if type(v) == dict or type(v) == ps.ParameterSet:
                    y[j] = (k, sort_deep_dict(v))
    return x


def which(program):
    '''
    taken from http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    '''
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

