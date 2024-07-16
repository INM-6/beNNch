"""
Copyright (C) 2023 Forschungszentrum Juelich GmbH, INM-6

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

import sys
import os
import shutil

from nbformat import v4 as nbf
import nbformat
import yaml
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.transforms as mtransforms

import bennchplot as bp


def plot_nodes(timer_file, save_path, plotting_style, detailed_timers):
    
    # If Docs is chosen as plotting style there are no detailed timers in the plot
    if plotting_style == 'Docs':
        detailed_timers = 'OFF'

    B = bp.Plot(data_file=timer_file, x_axis=['num_nodes'], time_scaling=1e3)

    if detailed_timers == 'ON':
        # Plotting with detailed timers
        widths = [1, 1]
        heights = [3, 1]
        fig = plt.figure(figsize=(12, 6), constrained_layout=True)
        spec = gridspec.GridSpec(ncols=2, nrows=2, figure=fig, width_ratios=widths, height_ratios=heights)

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
                            'wall_time_create+wall_time_connect',
                            'wall_time_sim', ],
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

        # Set the ylim for benchmark plots to make different runs of same model visually comparable
        if benchmark_model == "microcircuit":
            ax1.set_ylim(0, 15)
            ax2.set_ylim(0, 1.5)
        elif benchmark_model == "hpc_benchmark":
            ax1.set_ylim(0, 45)
            ax2.set_ylim(0, 45)
        elif benchmark_model == "multi_area":
            ax1.set_ylim(0, 1500)
            ax2.set_ylim(0, 150)
        else:
            pass # Let the ylim creation automatic for potential other models

        ax3.set_ylim(0, 100)

        for ax in [ax1, ax2, ax3]:
            ax.margins(x=0)
        for ax in [ax1, ax2]:
            B.simple_axis(ax)

    elif detailed_timers == 'OFF':
        if plotting_style == 'Standard':
        # Plotting without detailed timers
            fig = plt.figure(figsize=(12, 6), constrained_layout=True)
            spec = gridspec.GridSpec(ncols=1, nrows=1, figure=fig)

        elif plotting_style == 'Docs':
            fig = plt.figure(figsize=(12, 6), constrained_layout=False)
            spec = gridspec.GridSpec(ncols=2, nrows=1, figure=fig, hspace=0.2)

        ax1 = fig.add_subplot(spec[0, 0])

        trans = mtransforms.ScaledTranslation(-20 /
                                                72, 7 / 72, fig.dpi_scale_trans)

        B.plot_fractions(axis=ax1,
                        fill_variables=[
                            'wall_time_create+wall_time_connect',
                            'wall_time_sim', ],
                        interpolate=True,
                        step=None,
                        error=True)
        
        ax1.set_xlabel('Number of Nodes')
        ax1.set_ylabel(r'$T_{\mathrm{wall}}$ [s] for $T_{\mathrm{model}} =$'
                        + f'{np.unique(B.df.model_time_sim.values)[0]} s')

        handles1, labels1 = ax1.get_legend_handles_labels()
        ax1.legend(handles1[::-1], labels1[::-1])

        # Set the ylim for benchmark plots to make different runs of same model visually comparable
        if benchmark_model == "microcircuit":
            ax1.set_ylim(0, 15)
        elif benchmark_model == "hpc_benchmark":
            ax1.set_ylim(0, 45)
        elif benchmark_model == "multi_area":
            ax1.set_ylim(0, 1500)
        else:
            pass # Let the ylim creation automatic for potential other models
            
        ax1.margins(x=0)
        B.simple_axis(ax1)

        if plotting_style == "Docs":
            ax2 = fig.add_subplot(spec[0, 1])
            
            B.plot_main(quantities=['sim_factor'], axis=ax2,
                        error=True, fmt='-')

            ax2.set_xlabel('Number of Nodes')
            ax2.set_ylabel(r'real-time factor $T_{\mathrm{wall}}/$'
                            r'$T_{\mathrm{model}}$')

            # Set the xlim
            ax2.set_xlim(left=0)

            # Set the ylim for benchmark plots to make different runs of same model visually comparable
            if benchmark_model == "microcircuit":
                ax2.set_ylim(0, 1.)
            elif benchmark_model == "hpc_benchmark":
                ax2.set_ylim(0, 35)
            elif benchmark_model == "multi_area":
                ax2.set_ylim(0, 120)
            else:
                pass # Let the ylim creation automatic for potential other models

            ax2.hlines(y=1, xmin=0, xmax=np.max(B.df[B.x_axis].values), color='gray', linestyle='--')

            handles2, labels2 = ax2.get_legend_handles_labels()
            ax2.legend(handles2[::-1], labels2[::-1])
            
            ax2.margins(x=0)
            B.simple_axis(ax2)
            
            """
            ax1.text(0.0, 1.0, 'A', transform=ax1.transAxes + trans,
                    fontsize='medium', va='bottom', fontweight='bold')

            ax2.text(0.0, 1.0, 'B', transform=ax1.transAxes + trans,
                    fontsize='medium', va='bottom', fontweight='bold')
            """

    else:
        raise NameError('Detailed timers must either be set to ON or OFF')

    plt.savefig(save_path, dpi=400)


def plot_threads(timer_file, save_path):
    B = bp.Plot(data_file=timer_file, x_axis=['num_nvp'], time_scaling=1e3)

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
                     ])

    ax1.set_ylabel(r'$T_{\mathrm{wall}}$ [s] for $T_{\mathrm{model}} =$'
                    + f'{np.unique(B.df.model_time_sim.values)[0]} s')
    ax1.set_xlabel('number of vps')
    ax2.set_ylabel(r'$T_{\mathrm{wall}}$ [%]')
    B.merge_legends(ax1, ax2)

    plt.savefig(save_path, dpi=600)


display_plot_fn = r"""from IPython.display import Image
from IPython.display import HTML
def display_plot(result_table, plot_path, attributes, header, page_number=1):
    display_list = '<left><ul>\n'
    for attribute in attributes:
        value = os.popen(f'git annex metadata result_table --get {attribute}').read().strip()
        display_list += f'  <li>{attribute}: {value}</li>\n'
    display_list += '</ul></left>'

    display(HTML(f'<center><header>benchmark: {header}</header></center>'))
    display(Image(filename=plot_path))
    display(HTML(display_list))
    display(HTML(f'<center>page {page_number}</center>'))
"""

def make_notebook(plot_path: str, output_dir: str, notebook_filename: str, header: str, result_table, attributes_to_display):
    nb = nbf.new_notebook()
    cells = []

    codes = {
        'skip': [display_plot_fn],
        'slide': [f"display_plot('{result_table}', '{plot_path}', {attributes_to_display}, '{header}')"],
        'subslide': []
    }

    # page_number = 2
    # for timer_hash in timer_hashes[1:]:
    #     codes['subslide'].append(f"display_plot('{result_table}', '{plot_path}', {attributes_to_display}, '{header}', {page_number})")
    #     page_number += 1
    for key, code_list in codes.items():
        for code in code_list:
            cells.append(nbf.new_code_cell(code, metadata={"slideshow": {"slide_type": key}}))

    nb['cells'] = cells

    notebook_file = os.path.join(output_dir, notebook_filename)
    with open(notebook_file, 'w') as _:
        nbformat.write(nb, _)

    os.system(f'jupyter nbconvert --inplace --execute "{notebook_file}"')
    os.system(f'jupyter nbconvert --to slides "{notebook_file}" --output-dir="{output_dir}" '
              + '--TemplateExporter.exclude_input=True --SlidesExporter.reveal_transition="none"')
    os.remove(notebook_file)


if __name__ == '__main__':
    plotting_style = sys.argv[1]
    benchmark_model = sys.argv[2]
    detailed_timers = sys.argv[3]
    jube_results_dir = sys.argv[4]
    output_dir = sys.argv[5]
    pipeline_id = sys.argv[6]
    attributes_to_display = sys.argv[7:]

    plots_dir = os.path.join(output_dir, "plots")
    os.makedirs(plots_dir)

    for benchmark_id in os.listdir(jube_results_dir): 
        result_table = os.path.join(jube_results_dir, benchmark_id, "result", "result_table.dat")
        
        benchmark_name = os.popen(f"jube info {jube_results_dir} --id {benchmark_id}").read().split("\n")[1].split(" id")[0][2:]

        plot_nodes(timer_file=result_table, save_path=os.path.join(plots_dir, "nodes.png"), plotting_style, detailed_timers)
        make_notebook(os.path.join(plots_dir, "nodes.png"), output_dir, f"{pipeline_id}_{benchmark_name}_nodes.ipynb", f"{pipeline_id} - {benchmark_name} (nodes)", result_table, attributes_to_display)
        plot_threads(timer_file=result_table, save_path=os.path.join(plots_dir, "threads.png"))
        make_notebook(os.path.join(plots_dir, "threads.png"), output_dir, f"{pipeline_id}_{benchmark_name}_threads.ipynb", f"{pipeline_id} - {benchmark_name} (threads)", result_table, attributes_to_display)

