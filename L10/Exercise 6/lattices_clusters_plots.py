import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy.optimize import curve_fit

folder = 'res_autocorrelation/prova3.5/'
Temperatures = ['_0.5_T_c', '_T_c', '_2_T_c']
# 'T_c'
# types = ['energy', 'final_lattice_2D_Ising', 'magnetization']
# '.txt'

# make a color map of fixed colors
# cmap = colors.LinearSegmentedColormap.from_list('my_colormap', ['darkblue', 'darkred'], 256)
cmap = colors.LinearSegmentedColormap.from_list('my_colormap', ['blue', 'red'], 256)

# Lattice plots
LL = np.array([20, 30, 50, 100, 150, 200])

for k, j in enumerate(LL):
    for T in [0, 1_000_000]:
        file_name = folder+str(j)+'_'+str(T)+'_T_c_lattice.txt'
        data = np.loadtxt(file_name)
        plt.figure(T)
        im = plt.pcolormesh(data, edgecolors='k', cmap=cmap, vmin=-1, vmax=1)
        plt.colorbar(im, cmap=cmap)
        plt.axis('off')
        name = f'L = {j} - ' +r'$T_c$' + f' - Iteration = {T}'
        plt.title(name)
        ax = plt.gca()
        ax.set_aspect('equal')
        # plt.savefig(f'img/lattice{j}_iteration_{T}_lattice.png', dpi=300)
        plt.show(block=False)
    plt.show()
    
    # magn and energy plot:
    file_name = folder+str(j)+'_energy.txt'
    data = np.loadtxt(file_name)
    plt.plot(data[-1000:], label=r'$\mathcal{H}$', color='tab:blue')
    plt.hlines(np.mean(data), 0, 1000, ls='--', color='black')
    file_name = folder+str(j)+'_magnetization.txt'
    data = np.abs(np.loadtxt(file_name))
    plt.plot(data[-1000:], label=r'$\mathcal{m}$', color='tab:red')
    plt.hlines(np.mean(data), 0, 1000, ls='--', color='black')
    plt.legend()
    plt.title('Magnetization & Energy time series')
    plt.show()

    file_name = folder+str(j)+'_cluster_size.txt'
    data = np.loadtxt(file_name)
    plt.hist(data, label="data", density=True, bins=22, facecolor='mediumturquoise', edgecolor='black')
    plt.title(r'$\mathcal{S}(\beta)$')
    print(f'Average cluster size = {np.mean(data)}, L^2 = {j**2}')
    plt.show()
