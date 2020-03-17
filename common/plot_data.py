import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import os



def plot_data(df_list, save, title, axes, axes_label, log_scale, legend):
    '''
    Plot data

    Parameters
    ----------

    df_list : list
        list containing data

    save : string
        path or name to location where to save figure

    title : string

    axes : list
        what to plot as list of strings to fetch from dataframe

    axes_label : list
        label for axes as list strings

    log_scale : bool

    legend : list
        legned of plots as list of string
    '''


    fig, ax1 = plt.subplots(1, 1, figsize=(16, 9), dpi=100)

    plt.title(title)
    plt.legend(prop={'size': 9})

    ticks = np.int_([1, 2, 4, 8,  16, 32, 64, 128, 256])
    tick_labels = np.int_(ticks / 1)

    yticks = np.int_([1, 10, 100])

    plt.axvline(128, color='black')
    plt.axvline(64, color='grey')
    plt.axhline(1., color='r')

    yt_locs, yt_labels = plt.yticks()
    yt_labels = np.concatenate([yt_labels, [yt_labels[0]]])

    for i, df in enumerate(df_list):
        plt.plot(df[axes[0]], df[axes[1]], marker ='o', label = legend[i])
    y1 = 50
    M = 2 * 8
    y2 = y1 / M
    y = np.array([y1, y2])
    x = np.array([1,M])
    ax1.plot(x,y)


    x_label = axes_label[0]
    y_label = axes_label[1]

    ax1.set_xlabel(x_label)
    ax1.set_ylabel(y_label)

    if log_scale:
        ax1.set_xscale('log')
        ax1.set_yscale('log')

    ax1.set_xticks(ticks)
    ax1.set_xticklabels(tick_labels)
    ax1.minorticks_off()

    ax1.set_yticks(yticks)
    ax1.set_yticklabels(yticks)
    ax1.legend()

    plt.savefig(save + '.pdf')
    plt.show()



# get data and do some preprocessing (calculate real-time factor etc.)
basepath = '/home/kurth/work/repositories/projects/microcircuit_jube_benches/data/hamstein_diffs'
suffix ='result/internal_results_csv.dat'
data_ids = ['000000', '000001']

data_paths = [os.path.join(basepath, _id, suffix) for _id in data_ids]

data_frames = [pd.read_csv(path, delimiter = ';') for path in data_paths]

for df in data_frames:
    df['rtf'] = df['simulation_time_avg']/(df['Tsim']/1000.)
    df['virtual_procs'] = df['cpuspertask']*df['nodes']*df['ntasks']



# arguments for plotting routing
arguments = {'df_list' : data_frames,
             'save' : 'hambach45_46_mpich',
             'title' : 'nest2-14 with mpich on one node',
             'axes' : ['virtual_procs', 'rtf'],
             'axes_label' : ['Threads', r'Real time factor $T_{sim}/T_{bio}$'],
             'log_scale' : True,
             'legend' : ['hambach46', 'hambach45']}

plot_data(**arguments)
