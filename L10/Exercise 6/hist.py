import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

folder = 'MMC/'
T_c = 2./(np.log(1+np.sqrt(2))); 
T_vals = [T_c/2., 0.8*T_c, T_c, 1.5*T_c, 2.*T_c]
T_list = ['0.5_T_c_', '0.8_T_c_', 'T_c_', '1.5_T_c_', '2_T_c_']
T_latex = [r'$T_1 = 0.5 \cdot T_c$', r'$T_2 = 0.8 \cdot T_c$', 
           r'$T_3 = T_c$', r'$T_4 = 1.5 \cdot T_c$', r'$T_5 = 2\cdot T_c$']
keys = ['lattice', 'magnetization', 'energy', 'swapping_rates']
ex = '.txt'
L = 50

# energy distribution
n_bins = 50
colors = ['mediumseagreen', 'mediumturquoise', 'mediumblue', 
          'mediumpurple', 'mediumorchid', 'mediumslateblue']
colors = ['#067bc2', '#84bcda', '#ecc30b', '#f37748', '#d56062']
for k, T in enumerate(T_list):
    file_name = folder+T+keys[2]+ex
    data = np.loadtxt(file_name)
    plt.hist(data, label=f'{T_latex[k]}', density=True, 
             bins=n_bins, facecolor=colors[k], 
             edgecolor='black', alpha=0.7, lw=1.2)

plt.title(f'Energy distribution of adjacent chains')
plt.legend()
plt.grid(ls='--')
plt.xlabel(r'$\mathcal{H}$')
plt.ylabel(r'$\mathcal{P}\:(\mathcal{H})$')
plt.tight_layout()
plt.show()
