###########################
# Single Histogram Method #
###########################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

plt.rcParams.update({'font.size': 13})  
folder = 'res/MMC/50/'
T_c = 2./(np.log(1+np.sqrt(2)))
eps = 0.01
T_vals = np.array([T_c-2.*eps, T_c-eps, T_c, T_c+eps, T_c+2.*eps])
T_list = ["T_c_2eps_", "T_c_eps_", "T_c_", "T_c_p_eps_", "T_c_p_2eps_"]
T_latex = [r'$T_c - 2\varepsilon$', r'$T_c-\varepsilon$', 
           r'$T_c$', r'$T_c + \varepsilon$', 
           r'$T_c + 2\varepsilon$']

not_good_estimate = False

if not_good_estimate:
    folder = "../../L10/Exercise 4/MMC_no_swap/"
    T_c = 2./np.log(1+np.sqrt(2))
    T_vals = np.array([0.5*T_c, 0.8*T_c, T_c, 1.5*T_c, 2.*T_c])
    T_list = ["0.5_T_c_", "0.8_T_c_", "T_c_", "1.5_T_c_", "2_T_c_"]
    T_latex = [r'$0.5\cdot T_c$', r'$0.8\cdot T_c$', 
            r'$T_c$', r'$1.5\cdot T_c $', 
            r'$2\cdot T_c $']


keys = ['lattice', 'magnetization', 'energy', 'swapping_rates']
ex = '.txt'
colors = np.array(["#ddfff7", "#93e1d8", "#ffa69e", "#aa4465", "#460b2e"])
colors = colors[::-1]

# file_n ='res/shm/direct_average.txt'
# average_energy = np.loadtxt(file_n)

 
file_n ='res/shm.txt'
betas = 1./T_vals
L = 50
N = L**2

betas_range = np.linspace(min(betas), max(betas), num=50)
# betas_range = np.loadtxt('res/shm/betas.txt')

single_histo_method = []

verify_overlap = True
average_energy = []
if verify_overlap:
    for k, T in enumerate(T_list):
        data = np.loadtxt(f'{folder}{T}{keys[2]}{ex}')*N
        plt.hist(data, bins=50, color=colors[k], 
                 density=True, alpha=0.6, label=T_latex[k])
        average_energy.append(np.mean(data))
    plt.xlabel(r'$E_i$')
    plt.ylabel(r'$\rho(E_i)$')
    plt.grid(ls='--', alpha=0.5)
    plt.legend(frameon=False)
    plt.tight_layout()
    prefix = 'no_' if not_good_estimate else ''
    plt.savefig(f'img/{prefix}overlap.png', dpi=300)
    plt.show()

average_energy = np.array(average_energy)

for k, T in enumerate(T_vals):
    
    energy_estimates = np.loadtxt(f'res/shm/U_beta_{T_list[k]}{ex}')
    
    single_histo_method.append(energy_estimates)


for i in range(len(colors)):
    # label = r'' if i == 0 else ''
    plt.scatter(1./betas_range, single_histo_method[i]/N, s=60+i*10, 
                color=colors[i], marker='o', alpha=0.7,
                label=r'$U(T_'+f'{i+1}'+r')$')
    plt.scatter(T_vals, average_energy/N, s=100, 
                color='tab:red', marker='*')
    # cols = colors[colors != colors[i]]
    # for j in range(len(colors)-1):
    #     labs = np.array(range(1, len(colors)+1))
    #     labs = labs[labs != i+1]
    #     plt.scatter(i+1, single_histo_method[i][j], s=80, color=cols[j], marker='.', 
    #                 alpha=0.5) #, label=r'SHM @ $\beta_'+f'{labs[j]}'+r'$')
    
    plt.xlabel(r'$T$')
    plt.xticks(T_vals, T_latex)
    plt.ylabel(r'$U(T)$')
    plt.grid(ls='--', alpha=0.7)
    plt.legend()
    plt.tight_layout()
prefix = 'not_' if not_good_estimate else ''
plt.savefig(f'img/{prefix}good_estimate.png', dpi=300)
plt.show()
