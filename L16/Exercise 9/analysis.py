#######################
# Harmonic oscillator #
#######################

import numpy as np
import matplotlib.pyplot as plt
import os

skip = True

plt.rcParams.update({'font.size': 13})  
folder = 'res/'
T = 10
ext = '.dat'
types = ['Exact', 'Euler', 'Symp']
colors = ['lightseagreen', 'dodgerblue', 'navy']

if not skip:
    for dt in os.listdir(folder)[1:]:
        fig, ax = plt.subplots(1, 3, figsize=(12, 3))
        # plt.suptitle(r'$\Delta t =$'+f'{dt}')
        for k, kind in enumerate(types):
            data_x = np.loadtxt(f'{folder}{dt}/{kind}_x{ext}')
            data_p = np.loadtxt(f'{folder}{dt}/{kind}_p{ext}')
            ax[k].plot(data_x, data_p, label=kind, color=colors[k])
            ax[k].legend()
            ax[k].grid(ls='--', alpha=0.7)
            ax[k].set_xlabel(r'$x$')
            ax[k].set_ylabel(r'$p$')
        plt.tight_layout()
        plt.savefig(f'img/{dt}.png', dpi=300)
        plt.show()

# for dt in os.listdir(folder):
#     fig, ax = plt.subplots(1, 2, figsize=(12, 3))
#     # plt.suptitle(r'$\Delta t =$'+f'{dt}')
#     data_true_x = np.loadtxt(f'{folder}{dt}/{types[0]}_x{ext}')
#     data_true_p = np.loadtxt(f'{folder}{dt}/{types[0]}_p{ext}')
        
#     for k, kind in enumerate(types[1:]):
#         data_x = np.loadtxt(f'{folder}{dt}/{kind}_x{ext}')
#         data_p = np.loadtxt(f'{folder}{dt}/{kind}_p{ext}')
#         ax[0].plot(data_true_x, 'r-', label='Exact')
#         ax[0].plot(data_x, label=kind, color='b')
#         ax[1].plot(data_true_p, 'r-', label='Exact')
#         ax[1].plot(data_p, label=kind, color='b')
#         plt.show()



# Residual plot to evaluate error
skip = True
if not skip:
    for j, dt in enumerate(os.listdir(folder)[1:]):
        loc = 'upper left' # if j == 0 else 'upper right'
        fig, ax = plt.subplots(2, 2, figsize=(12, 6))
        # plt.suptitle(r'$\Delta t =$'+f'{dt}')
        data_true_x = np.loadtxt(f'{folder}{dt}/{types[0]}_x{ext}')
        data_true_p = np.loadtxt(f'{folder}{dt}/{types[0]}_p{ext}')
        
        # remove_x = data_true_x != 0
        # remove_p = data_true_p != 0
        max_x, max_p, min_x, min_p = np.zeros(4)

        for k, kind in enumerate(types[1:]):
            data_x = np.loadtxt(f'{folder}{dt}/{kind}_x{ext}')
            data_p = np.loadtxt(f'{folder}{dt}/{kind}_p{ext}')
            
            error_x = data_x-data_true_x
            error_p = data_p-data_true_p
            

            ax[k, 0].plot(np.arange(0, T+float(dt), float(dt)), error_x, label=kind, color=colors[k])
            ax[k, 0].legend(loc=loc)
            ax[k, 0].grid(ls='--', alpha=0.7)
            ax[k, 0].set_xlabel(r'$t$')
            ax[k, 0].set_ylabel(r'$\varepsilon_x$')
            # chi = round(sum((data_x-data_true_x)**2)/len(data_x), 4)
            # chi = chi if len(str(chi)) == 6 else str(chi).zfill(6)[::-1]
            # ax[k, 0].text(3, -0.0025 if j == 0 else -0.025, r'$\chi^2=$'+f'{chi}')

            min_x = min_x if min_x < ax[k, 0].get_ylim()[0] else ax[k, 0].get_ylim()[0]
            max_x = max_x if max_x > ax[k, 0].get_ylim()[1] else ax[k, 0].get_ylim()[1]
            
            ax[k, 1].plot(np.arange(0, T+float(dt), float(dt)), error_p, label=kind, color=colors[k])
            ax[k, 1].legend(loc=loc)
            ax[k, 1].grid(ls='--', alpha=0.7)
            ax[k, 1].set_xlabel(r'$t$')
            ax[k, 1].set_ylabel(r'$\varepsilon_p$')
            # chi = round(sum((data_p-data_true_p)**2)/len(data_x), 4)
            # chi = chi if len(str(chi)) == 6 else str(chi).zfill(6)[::-1]
            # ax[k, 1].text(1, -0.0025 if j == 0 else -0.025, r'$\chi^2=$'+f'{chi}')

            min_p = min_p if min_p < ax[k, 1].get_ylim()[0] else ax[k, 1].get_ylim()[0]
            max_p = max_p if max_p > ax[k, 1].get_ylim()[1] else ax[k, 1].get_ylim()[1]

        for i, axs in enumerate(ax.flat):
            if i % 2 == 0:
                axs.set_ylim((min_x, max_x))
            else:
                axs.set_ylim((min_p, max_p))
        plt.tight_layout()
        plt.savefig(f'img/error_{dt}.png', dpi=300)
        plt.show()

# Hamiltonian plot to evaluate error
def Hamiltonian(x, p, dt=0):
    H = (p**2 + x**2)/2
    if dt:
        return H - x*p*dt/2 
    return H

from matplotlib.ticker import MaxNLocator
skip = True
if not skip:   
    for j, dt in enumerate(os.listdir(folder)[1:]):
        fig, ax = plt.subplots(1, 2, figsize=(12, 4))
        # plt.suptitle(r'$\Delta t =$'+f'{dt}')
        data_true_x = np.loadtxt(f'{folder}{dt}/{types[0]}_x{ext}')
        data_true_p = np.loadtxt(f'{folder}{dt}/{types[0]}_p{ext}')
        
        H_exact = Hamiltonian(data_true_x, data_true_p)

        max_y, min_y = 0, 100

        for k, kind in enumerate(types[1:]):
            data_x = np.loadtxt(f'{folder}{dt}/{kind}_x{ext}')
            data_p = np.loadtxt(f'{folder}{dt}/{kind}_p{ext}')
            
            H_first = Hamiltonian(data_x, data_p)

            ax[k].plot(np.arange(0, T+float(dt), float(dt)), H_exact, 
                       label='Exact '+r'$\mathcal{H}$', 
                       color=colors[0])
            ax[k].plot(np.arange(0, T+float(dt), float(dt)), H_first, 
                       label=kind+r' $\mathcal{H}$', color=colors[1])

            y_val = 0.5025 if j == 0 else 0.525
            # ax[k].text(1.5, y_val, r'$\chi^2_{\mathcal{H}}=$'+f'{round(np.sum((H_exact-H_first)**2)/len(H_exact), 2)}')
            
            if kind == 'Symp':
                H_shadow = Hamiltonian(data_x, data_p, float(dt))
                ax[k].plot(np.arange(0, T+float(dt), float(dt)), H_shadow,
                           label=f'Shadow '+ r'$\mathcal{H}$'+f"'", color=colors[2])
                y_val = 0.5020 if j == 0 else 0.520
                # ax[k].text(1.5, y_val, r"$\chi^2_{\mathcal{H}'}=$"+f'{round(np.sum((H_exact-H_shadow)**2)/len(H_exact), 2)}')
            
            ax[k].legend(loc='upper left')
            ax[k].grid(ls='--', alpha=0.7)
            ax[k].set_xlabel(r'$t$')
            ax[k].set_ylabel(r'$\mathcal{H}$')

            min_y = min_y if min_y < ax[k].get_ylim()[0] else ax[k].get_ylim()[0]
            max_y = max_y if max_y > ax[k].get_ylim()[1] else ax[k].get_ylim()[1]
            
            
        for i, axs in enumerate(ax.flat):
            axs.set_ylim((min_y, max_y))
            axs.yaxis.set_major_locator(MaxNLocator(nbins=4))

            
        plt.tight_layout()
        plt.savefig(f'img/Hamiltonian_{dt}.png', dpi=300)
        plt.show()

############################
# Velocity Verlet analysis #
############################
        
types = ['Exact_omega', 'Velocity_Verlet', 'Beeman']
labels = ['Exact', 'Velocity\nVerlet', 'Beeman']
dt = '0.001'

skip = True
if not skip:
    fig, ax = plt.subplots(1, 3, figsize=(12, 3))
    for k, kind in enumerate(types):
        data_x = np.loadtxt(f'{folder}{dt}/{kind}_x{ext}')
        data_p = np.loadtxt(f'{folder}{dt}/{kind}_p{ext}')
        ax[k].plot(data_x, data_p, label=labels[k], color=colors[k])
        ax[k].legend()
        ax[k].grid(ls='--', alpha=0.7)
        ax[k].set_xlabel(r'$x$')
        ax[k].set_ylabel(r'$p$')
    plt.tight_layout()
    plt.savefig(f'img/VV.png', dpi=300)
    plt.show()

# Residual plot to evaluate error
skip = True
if not skip:
    loc = 'upper left' # if j == 0 else 'upper right'
    fig, ax = plt.subplots(2, 2, figsize=(12, 6))
    data_true_x = np.loadtxt(f'{folder}{dt}/{types[0]}_x{ext}')
    data_true_p = np.loadtxt(f'{folder}{dt}/{types[0]}_p{ext}')
    
    max_x, max_p, min_x, min_p = np.zeros(4)

    for k, kind in enumerate(types[1:]):
        data_x = np.loadtxt(f'{folder}{dt}/{kind}_x{ext}')
        data_p = np.loadtxt(f'{folder}{dt}/{kind}_p{ext}')
        
        error_x = data_x-data_true_x
        error_p = data_p-data_true_p

        ax[k, 0].plot(np.arange(0, T+float(dt), float(dt)), error_x, label=labels[k+1], color=colors[k])
        ax[k, 0].legend(loc=loc)
        ax[k, 0].grid(ls='--', alpha=0.7)
        ax[k, 0].set_xlabel(r'$t$')
        ax[k, 0].set_ylabel(r'$\varepsilon_x$')
        # ax[k, 0].text(2 , -2e-7, r'$\chi^2=$'+f'{round(sum(error_x**2)/len(data_x), 4)}')

        min_x = min_x if min_x < ax[k, 0].get_ylim()[0] else ax[k, 0].get_ylim()[0]
        max_x = max_x if max_x > ax[k, 0].get_ylim()[1] else ax[k, 0].get_ylim()[1]
        
        ax[k, 1].plot(np.arange(0, T+float(dt), float(dt)), error_p, label=labels[k+1], color=colors[k])
        ax[k, 1].legend(loc=loc)
        ax[k, 1].grid(ls='--', alpha=0.7)
        ax[k, 1].set_xlabel(r'$t$')
        ax[k, 1].set_ylabel(r'$\varepsilon_p$')
        # chi = chi if len(str(chi)) == 6 else str(chi).zfill(6)[::-1]
        # ax[k, 1].text(2, -2e-7, r'$\chi^2=$'+f'{round(sum(error_p**2)/len(data_x), 4)}')

        min_p = min_p if min_p < ax[k, 1].get_ylim()[0] else ax[k, 1].get_ylim()[0]
        max_p = max_p if max_p > ax[k, 1].get_ylim()[1] else ax[k, 1].get_ylim()[1]

    for i, axs in enumerate(ax.flat):
        if i % 2 == 0:
            axs.set_ylim((min_x, max_x))
        else:
            axs.set_ylim((min_p, max_p))
    plt.tight_layout()
    plt.savefig(f'img/VV_error_{dt}.png', dpi=300)
    plt.show()

# Energy plot to evaluate error
skip = True
omega = 1.
if not skip:   
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    # plt.suptitle(r'$\Delta t =$'+f'{dt}')
    data_true_x = np.loadtxt(f'{folder}{dt}/{types[0]}_x{ext}')
    data_true_p = np.loadtxt(f'{folder}{dt}/{types[0]}_p{ext}')
    
    H_exact = ((data_true_x*omega)**2 + (data_true_p)**2)/2.
    H_exact_rel = (H_exact - H_exact[0])/H_exact[0]
    max_y, min_y = 0, 100

    for k, kind in enumerate(types[1:]):
        data_x = np.loadtxt(f'{folder}{dt}/{kind}_x{ext}')
        data_p = np.loadtxt(f'{folder}{dt}/{kind}_p{ext}')
        
        H_kind = ((data_x*omega)**2 + (data_p)**2)/2.
        H_kind_rel = (H_kind-H_kind[0])/H_kind[0]

        ax[k].plot(np.arange(0, T+float(dt), float(dt)), H_exact_rel, 
                    label='Exact '+r'$\Delta\mathcal{H}$', lw = 3,
                    color=colors[0])
        ax[k].plot(np.arange(0, T+float(dt), float(dt)), H_kind_rel, 
                    label=labels[k+1]+r' $\Delta\mathcal{H}$', color=colors[k+1])

        y_val = 2e-8 if k == 0 else -1e-7
        # ax[k].text(2, y_val, r'$\chi^2_{\mathcal{H}}=$'+f'{round(np.sum((H_exact-H_kind)**2)/len(H_exact), 4)}')
        
        ax[k].legend(loc='upper right' if k == 0 else 'lower right')
        ax[k].grid(ls='--', alpha=0.7)
        ax[k].set_xlabel(r'$t$')
        ax[k].set_ylabel(r'$\Delta\mathcal{H}$')

        min_y = min_y if min_y < ax[k].get_ylim()[0] else ax[k].get_ylim()[0]
        max_y = max_y if max_y > ax[k].get_ylim()[1] else ax[k].get_ylim()[1]
        
        
    for i, axs in enumerate(ax.flat):
        axs.set_ylim((min_y, max_y))
        axs.yaxis.set_major_locator(MaxNLocator(nbins=4))

        
    plt.tight_layout()
    plt.savefig(f'img/VV_Hamiltonian_{dt}.png', dpi=300)
    plt.show()


