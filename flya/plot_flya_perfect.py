import numpy as  np
import matplotlib as mpl
mpl.use("Agg")
import astropy.table as tab
import matplotlib.pyplot as plt

def plot_flya(simname):
    # setting the figure
    font = {'family': 'serif', 'weight': 'normal', 'size': 13}
    plt.rc('font', **font)
    mpl.rcParams['axes.linewidth'] = 1.5

    out_fig_name = 'tau_perfect_{}.pdf'.format(simname)
    figure_size = [7, 7]
    fig, ax = plt.subplots(1, 1, figsize=(figure_size[0], figure_size[1]))

    path = '/home/vikram/flya/flya/data'
    tau_list = list([0, 0.01, 0.05, 0.1])
    c_list = ['red', 'green', 'blue', 'orange']

    for Gamma12, ls in zip([0.05, 0.075, 0.1], [':', '-', '--']):
        for tau_low, c in zip(tau_list, c_list):
            label = 'G{:0.3f}-tau{:0.2f}'.format(Gamma12, tau_low)
            file_name = path + '/' + 'dflya_{}_Gamma_{:0.3f}_taulim_{:0.2f}_x.fits'.format(simname, Gamma12, tau_low)
            data = tab.Table.read(file_name)
            ax.plot(data['tau'], data['frac'], label=label, linewidth=2, color=c, linestyle=ls)

    ax.set_ylabel(r'f$_{\rm Ly \alpha}$')
    ax.set_xlabel(r'$\tau_{\rm lim}$')

    ax.legend(loc='best', fontsize=10)
    ax.set_ylim(0.15, 0.55)
    # decorating the plot
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(1.5)

    ax.grid()

    fig.savefig(out_fig_name, bbox_inches='tight')


plot_flya('tng')
plot_flya('ill')