"""
beNNch - Unified execution, collection, analysis and
comparison of neural network simulation benchmarks.
Copyright (C) 2021 Forschungszentrum Juelich GmbH, INM-6

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with
this program. If not, see <https://www.gnu.org/licenses/>.

SPDX-License-Identifier: GPL-3.0-or-later
"""

import numpy as np
import bennchplot as bp
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.transforms as mtransforms
import tol_colors
import os


def plot(scaling_type, timer_hash, timer_file, save_path):

    if scaling_type == 'nodes':
        args = {
            'data_file': timer_file,
            'x_axis': ['num_nodes'],
            'time_scaling': 1e3,
            'detailed_timers': True,
        }

        # Instantiate class
        B = bp.Plot(**args)

        # Plotting
        widths = [1, 1]
        heights = [3, 1]
        fig = plt.figure(figsize=(12, 6), constrained_layout=True)
        spec = gridspec.GridSpec(ncols=2, nrows=2, figure=fig,
                                 width_ratios=widths,
                                 height_ratios=heights)

        ax1 = fig.add_subplot(spec[:, 0])
        ax2 = fig.add_subplot(spec[0, 1])
        ax3 = fig.add_subplot(spec[1, 1])

        trans = mtransforms.ScaledTranslation(-20 /
                                              72, 7 / 72, fig.dpi_scale_trans)
        ax1.text(0.0, 1.0, 'A', transform=ax1.transAxes + trans,
                 fontsize='medium', va='bottom', fontweight='bold')
        ax2.text(0.0, 1.0, 'B', transform=ax2.transAxes + trans,
                 fontsize='medium', va='bottom', fontweight='bold')
        ax3.text(0.0, 1.0, 'C', transform=ax3.transAxes + trans,
                 fontsize='medium', va='bottom', fontweight='bold')

        B.plot_fractions(axis=ax1,
                         fill_variables=[
                             'time_construction_create+time_construction_connect',
                             'time_simulate', ],
                         interpolate=True,
                         step=None,
                         error=True)
        B.plot_main(quantities=['sim_factor'], axis=ax2,
                    error=True)
        B.plot_fractions(axis=ax2,
                         fill_variables=[
                             'phase_update_factor',
                             'phase_collocate_factor',
                             'phase_communicate_factor',
                             'phase_deliver_factor'
                         ])
        B.plot_fractions(axis=ax3,
                         fill_variables=[
                             'frac_phase_update',
                             'frac_phase_collocate',
                             'frac_phase_communicate',
                             'frac_phase_deliver'
                         ])

        ax1.set_xlabel('Number of Nodes')
        ax1.set_ylabel(r'$T_{\mathrm{wall}}$ [s] for $T_{\mathrm{model}} =$'
                       + f'{np.unique(B.df.model_time_sim.values)[0]} s')
        ax2.set_ylabel(r'real-time factor $T_{\mathrm{wall}}/$'
                       r'$T_{\mathrm{model}}$')
        ax3.set_xlabel('Number of Nodes')
        ax3.set_ylabel(r'relative $T_{\mathrm{wall}}$ [%]')

        handles1, labels1 = ax1.get_legend_handles_labels()
        handles2, labels2 = ax2.get_legend_handles_labels()

        ax1.legend(handles1[::-1], labels1[::-1])
        ax2.legend(handles2[::-1], labels2[::-1], loc='upper right')

        ax3.set_ylim(0, 100)
        ax1.set_ylim(0, 65)
        ax2.set_ylim(0, 1.3)

        for ax in [ax1, ax2, ax3]:
            ax.margins(x=0)
        for ax in [ax1, ax2]:
            B.simple_axis(ax)

        plt.savefig(f'{save_path}/{timer_hash}.png', dpi=400)

    elif scaling_type == 'threads':

        args = {
            'data_file': timer_file,
            'x_axis': ['num_nvp'],
            'time_scaling': 1e3
        }

        # Instantiate class
        B = bp.Plot(**args)

        # Plotting
        widths = [1]
        heights = [3, 1]
        fig = plt.figure(figsize=(6, 6), constrained_layout=True)
        spec = gridspec.GridSpec(ncols=1, nrows=2, figure=fig,
                                 width_ratios=widths,
                                 height_ratios=heights)

        ax1 = fig.add_subplot(spec[0, :])
        ax2 = fig.add_subplot(spec[1, :])

        B.plot_main(quantities=['sim_factor'], axis=ax1, log=(False, True))
        B.plot_fractions(axis=ax2,
                         fill_variables=[
                             'frac_phase_update',
                             'frac_phase_collocate',
                             'frac_phase_communicate',
                             'frac_phase_deliver'
                         ],
                         )

        ax1.set_ylabel(r'$T_{\mathrm{wall}}$ [s] for $T_{\mathrm{model}} =$'
                       + f'{np.unique(B.df.model_time_sim.values)[0]} s')
        ax1.set_xlabel('number of vps')
        ax2.set_ylabel(r'$T_{\mathrm{wall}}$ [%]')
        B.merge_legends(ax1, ax2)

        plt.savefig(f'{save_path}/{timer_hash}.png', dpi=600)


def plot_comparison(scaling_type, timer_files, save_path, colors=None):

    if colors is None:
        vibrant = tol_colors.tol_cset('vibrant')
        colors = [vibrant.blue, vibrant.orange, vibrant.red, vibrant.teal, vibrant.magenta, vibrant.cyan]

    if len(timer_files) > len(colors):
        raise NotImplementedError('Number of colors < number of lines in plot. Please provide a longer list of colors.')

    fig, ax = plt.subplots(figsize=(6, 6), constrained_layout=True)

    if scaling_type == 'nodes':
        xaxis = ['num_nodes']
        xaxis_label = 'Number of Nodes'
    elif scaling_type == 'threads':
        xaxis = ['num_nvp']
        xaxis_label = 'Number of Virtual Processes'
    for i, timer_file in enumerate(timer_files):

        args = {
            'data_file': timer_file,
            'x_axis': xaxis,
            'time_scaling': 1e3,
            'detailed_timers': False
        }
        timer_hash = timer_file.split('/')[-1].split('.')[0]
        file_path = os.popen(f"find . -name '*{timer_hash}.csv'").read().strip()
        version = os.popen(f'git annex metadata {file_path} --get simulator-version').read().strip()

        # Instantiate class
        B = bp.Plot(**args)
        B.plot_main(quantities=['sim_factor'], axis=ax, error=True, label=version, color=colors[i])
        B.simple_axis(ax)
    ax.legend()
    ax.set_ylabel(r'real-time factor $T_{\mathrm{wall}}/T_{\mathrm{model}}$')
    ax.set_xlabel(xaxis_label)

    plt.savefig(f'{save_path}/comparison.png', dpi=600)
