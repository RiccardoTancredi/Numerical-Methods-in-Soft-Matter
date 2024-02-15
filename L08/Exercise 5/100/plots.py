import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
# import seaborn as sns

folder = 'res/'
Temperatures = ['0_1_', '0_9_', '2_']
# 'T_c'
types = ['energy', 'final_lattice_2D_Ising', 'magnetization']
# '.txt'

# make a color map of fixed colors
# cmap = colors.ListedColormap(['darkblue', 'darkred'])
cmap = colors.LinearSegmentedColormap.from_list('my_colormap', ['blue', 'red'], 256)
cmap2 = colors.LinearSegmentedColormap.from_list('my_colormap2', ['blue', 'red'], 256)

# Lattice plots
lattice_plots = True
L = 100
if lattice_plots:
    for T in Temperatures:
        file_name = folder+T+'T_c_'+types[1]+'.txt'
        data = np.loadtxt(file_name)
        plt.clf()
        im = plt.pcolormesh(data, edgecolors='k', cmap=cmap, vmin=-1, vmax=1)
        plt.colorbar(im, cmap=cmap)
        plt.axis('off')
        name = '.'.join(T.split("_"))[:-1] +r'$\cdot T_c$'
        plt.title(r'$L='+f'{L}'+r'$')
        ax = plt.gca()
        ax.set_aspect('equal')
        plt.savefig(f'img/lattice_{T}T_c.png', dpi=300)
        # plt.show()

# Energy & magnetization plots
colors = ["#4B88A2", "#BB0A21"]
labels = [r"$\mathcal{H}(\sigma)$", r"$m$"]
equilibrium_times = []
plt.clf()
for k, T in enumerate(Temperatures):
    file_name_e = folder+T+'T_c_'+types[0]+'.txt'
    file_name_m = folder+T+'T_c_'+types[2]+'.txt'
    energy = np.loadtxt(file_name_e)
    magn = np.loadtxt(file_name_m)
    # Find the equilibium time
    ind_e = np.where(abs(energy) > abs(energy.mean()))[0][0]
    ind_m = np.where(abs(magn) > abs(magn.mean()))[0][0]
    equilibrium_times.append(max([ind_e, ind_m]))
    # print(equilibrium_times)

    plt.plot(energy[0:equilibrium_times[k]+1000], c=colors[0], label=labels[0])
    plt.plot(magn[0:equilibrium_times[k]+1000], c=colors[1], label=labels[1])
    # if T != '2_':
    plt.vlines(equilibrium_times[k], energy[equilibrium_times[k]], 
            magn[equilibrium_times[k]], color='k', ls='--', label=r'$\tau_{eq}$')
        
    plt.xlabel('MC steps per lattice site')
    plt.ylabel('Observable')
    name = '.'.join(T.split("_"))[:-1] +r'$\cdot T_c$'
    plt.title(name)
    plt.legend()
    plt.grid(ls='--', alpha=0.5)
    plt.savefig(f'img/energy_magn_{T}T_c.png', dpi=300)
    # plt.show()
    plt.clf()

# All in one
maximum = max(equilibrium_times)
for k, T in enumerate(Temperatures):
    file_name_e = folder+T+'T_c_'+types[0]+'.txt'
    file_name_m = folder+T+'T_c_'+types[2]+'.txt'
    energy = np.loadtxt(file_name_e)
    magn = np.loadtxt(file_name_m)

    sub_lab1 = labels[0] if k == 0 else ""
    sub_lab2 = labels[1] if k == 0 else ""
    sub_lab3 = r'$\tau_{eq}$' if k == 0 else ""
    plt.plot(energy[0:maximum+1000], c=colors[0], label=sub_lab1)
    plt.plot(magn[0:maximum+1000], c=colors[1], label=sub_lab2)
    plt.vlines(equilibrium_times[k], energy[equilibrium_times[k]], 
            magn[equilibrium_times[k]], color='k', ls='--', label=sub_lab3)

    plt.xlabel('MC steps per lattice site')
    plt.ylabel('Observable')
    name = '.'.join(T.split("_"))[:-1] +r'$\cdot T_c$'
    plt.title(f'2D Ising model - L = {L}')
    plt.text(2000, -1.9, r'$T = 0.2\cdot T_c$')
    plt.text(8000, -1.55, r'$T = 0.9\cdot T_c$')
    plt.text(7000, 0.1, r'$T = 2\cdot T_c$')
    plt.legend(bbox_to_anchor=(0.3, 0.2))
    plt.grid(ls='--', alpha=0.5)
    # plt.show()
    plt.savefig(f'img/energy_magn', dpi=300)
    plt.clf()

    # Save T_eq for further analysis
    np.savetxt('t_eq.txt', equilibrium_times)

# many plots
folder += 'many/'
zeros = ['0_', '1000_', '10000_', '100000_', '1000000_']
for zer in zeros:
    file_name = folder+zer+'lattice_2D_Ising.txt'
    data = np.loadtxt(file_name)
    plt.clf()
    im = plt.pcolormesh(data, edgecolors='k', cmap=cmap, vmin=-1, vmax=1)
    plt.colorbar(im, cmap=cmap)
    plt.axis('off')
    n = zer.split("_")[0].count('0')
    name = r'$10^' f'{n}'+r'$' + f' MC step per lattice - L = {L}, ' + r'$T=0.9\cdot T_c$'
    plt.title(name)
    ax = plt.gca()
    ax.set_aspect('equal')
    plt.savefig(f'img/{zer}lattice.png', dpi=300)
    # plt.show()

plt.clf()
fig, ax = plt.subplots(1, len(zeros), figsize=(12, 3))
for i, zer in enumerate(zeros):
    file_name = folder+zer+'lattice_2D_Ising.txt'
    data = np.loadtxt(file_name)
    im = ax[i].pcolormesh(data, cmap=cmap2, vmin=-1, vmax=1)
    ax[i].axis('off')
    n = zer.split("_")[0].count('0')
    name = r'$10^' f'{n}'+r'$' # + f' MC step per lattice - L = {L}' + r'$0.9\cdot T_c$'
    ax[i].set_title(name)
    ax[i].set_aspect('equal')

plt.draw()
p0 = ax[0].get_position().get_points().flatten()
p1 = ax[1].get_position().get_points().flatten()
p2 = ax[2].get_position().get_points().flatten()
p3 = ax[3].get_position().get_points().flatten()
p4 = ax[4].get_position().get_points().flatten()

# ax_cbar = fig.add_axes([p0[0], 1, p4[2]-p0[0], 0.15])
# plt.colorbar(im, cax=ax_cbar, cmap=cmap2, orientation='horizontal')
# ax_cbar = fig.add_axes([p0[0], 1, p4[2]-p0[0], 0.05])
# plt.colorbar(im, cax=ax_cbar, cmap=cmap2, orientation='horizontal')
plt.savefig(f'img/many_lattices.png', dpi=100)
# plt.show()