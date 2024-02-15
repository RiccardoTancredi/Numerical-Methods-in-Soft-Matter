import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy.optimize import curve_fit

folder = 'step_by_step/'
Temperatures = ['_T_c', '_2_T_c', '_0.5_T_c']
# 'T_c'
# types = ['energy', 'final_lattice_2D_Ising', 'magnetization']
# '.txt'

# make a color map of fixed colors
# cmap = colors.ListedColormap(['darkblue', 'darkred'])
cmap = colors.LinearSegmentedColormap.from_list('my_colormap', ['blue', 'red'], 256)

# Lattice plots
lattice_plots = True

if lattice_plots:
    for T in range(0, 151, 10):
        file_name = folder+str(T)+'_T_c_lattice.txt'
        data = np.loadtxt(file_name)
        plt.figure(T)
        im = plt.pcolormesh(data, edgecolors='k', cmap=cmap, vmin=-1, vmax=1)
        plt.gca().invert_yaxis()
        plt.colorbar(im, cmap=cmap)
        plt.axis('off')
        # name = '.'.join(T.split("_"))[:-1] +r'$\cdot T_c$'
        plt.title(T)
        ax = plt.gca()
        ax.set_aspect('equal')
        # plt.savefig(f'img/lattice_{T}_lattice.png', dpi=300)
        plt.show(block=False)
plt.show()
