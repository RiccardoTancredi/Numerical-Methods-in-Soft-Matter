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

# Lattice plots
lattice_plots = True
L = 50
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
    # if T != '2_':
    plt.vlines(equilibrium_times[k], energy[equilibrium_times[k]], 
            magn[equilibrium_times[k]], color='k', ls='--', label=sub_lab3)

plt.xlabel('MC steps per lattice site')
plt.ylabel('Observable')
name = '.'.join(T.split("_"))[:-1] +r'$\cdot T_c$'
plt.title(f'2D Ising model - L = {L}')
plt.text(1000, -1.85, r'$T = 0.2\cdot T_c$')
plt.text(3000, -1.4, r'$T = 0.9\cdot T_c$')
plt.text(4000, 0.2, r'$T = 2\cdot T_c$')
plt.legend(bbox_to_anchor=(0.3, 0.2))
plt.grid(ls='--', alpha=0.5)
# plt.show()
plt.savefig(f'img/energy_magn', dpi=300)
plt.clf()

# Save T_eq for further analysis
np.savetxt('t_eq.txt', equilibrium_times)