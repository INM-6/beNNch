import benchplot as bp
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec


def plot(model, timer_hash, timer_path):

    if model in ['mam', 'MAM', 'multi-area-model']:
        args = {
            'data_hash': timer_hash,
            'data_path': timer_path,
            'x_axis': ['num_nodes'],
            'time_scaling': 1e3
        }

        # Instantiate class
        B = bp.BenchPlot(**args)

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

        B.plot_fractions(axis=ax1,
                         fill_variables=['wall_time_sim',
                                         'wall_time_create+wall_time_connect'],
                         interpolate=True,
                         step=None,
                         error=True)
        B.plot_main(quantities=['sim_factor', 'phase_total_factor'], axis=ax2,
                    error=True)
        B.plot_fractions(axis=ax3,
                         fill_variables=[
                             'frac_phase_communicate',
                             'frac_phase_update',
                             'frac_phase_deliver',
                             'frac_phase_collocate'
                         ])

        ax1.set_xlabel('Number of Nodes')
        ax1.set_ylabel('wall time [s]')
        ax2.set_ylabel(r'real-time factor $T_{\mathrm{wall}}/$'
                       r'$T_{\mathrm{model}}$')
        ax3.set_xlabel('Number of Nodes')

        ax1.legend()
        B.merge_legends(ax2, ax3)

        plt.savefig(f'{timer_path}/{timer_hash}.pdf')

    elif model in ['mc', 'MC', 'microcircuit']:

        args = {
            'data_hash': timer_hash,
            'data_path': timer_path,
            'x_axis': ['num_nvp'],
            'time_scaling': 1e3
        }

        # Instantiate class
        B = bp.BenchPlot(**args)

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
                             'frac_phase_communicate',
                             'frac_phase_update',
                             'frac_phase_deliver',
                             'frac_phase_collocate'
                         ],
                         )

        ax1.set_ylabel(r'real-time factor $T_{\mathrm{wall}}/$'
                       r'$T_{\mathrm{model}}$')
        ax1.set_xlabel('number of vps')
        ax1.legend()
        ax2.set_ylabel(r'relative wall time $[\%]$')
        B.merge_legends(ax1, ax2)

        plt.savefig(f'{timer_path}/{timer_hash}.pdf')
