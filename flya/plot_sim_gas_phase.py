import numpy as  np
import matplotlib as mpl
mpl.use("Agg")
import astropy.table as tab
import matplotlib.pyplot as plt

# setting the figure
font = {'family': 'serif', 'weight': 'normal', 'size': 13}
plt.rc('font', **font)
mpl.rcParams['axes.linewidth'] = 1.5

out_fig_name = 'gas_phase.pdf'
figure_size = [7, 7]
fig, ax = plt.subplots(1, 1, figsize=(figure_size[0], figure_size[1]))

path = '/home/vikram/flya/flya/data'
Teme_list = [5e4, 1e5, 5e5, 1e6, 5e6]
c_list = ['red', 'green', 'blue', 'orange', 'magenta']
for Temp, c in zip(Teme_list, c_list):

    simname = 'tng'
    file = path +'/' +'sim_gas_phase_{}_T_{:.0f}.fits'.format(simname, Temp)
    label = 'Tmax = {:0.0e} K'.format(Temp)
    data = tab.Table.read(file)
    ax.plot(data['oden'], data['diff'],  label=label, linewidth = 2, color = c)

    simname = 'ill'
    file = path +'/' +'sim_gas_phase_{}_T_{:.0f}.fits'.format(simname, Temp)
    label = ''
    data = tab.Table.read(file)
    ax.plot(data['oden'], data['diff'],  label=label, linewidth = 2, color = c, linestyle  = '--')

ax.set_ylabel(r'f$_{\rm Ly \alpha}$')
ax.set_xlabel(r'$\Delta$')

ax.legend( loc = 'best', fontsize = 10)
ax.set_ylim(0.15, 0.9)
# decorating the plot
for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(1.5)

ax.grid()

fig.savefig(out_fig_name, bbox_inches='tight')

