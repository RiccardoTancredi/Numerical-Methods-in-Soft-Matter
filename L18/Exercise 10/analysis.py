###############
# Thermostats #
###############

import numpy as np
import matplotlib.pyplot as plt
import os


plt.rcParams.update({'font.size': 13})  
folder = 'res/Velocity_Verlet/'
saveimg = '../img/'
sigma = 1.
sigma_cut = np.arange(1., 4.2, 0.2*sigma)
sigma_cut = [round(el, 2) for el in sigma_cut]

ext = '.dat'
types = ['energy', 'kinetic_energy', 'data', 'velocity']
colors = ['lightseagreen', 'dodgerblue', 'navy']

def compute_distances(data):
    return np.array([np.array([np.sqrt(np.sum((p-q - np.round((p-q)/10.)*10.)**2)) 
                               for j, q in enumerate(data)])
                               for i, p in enumerate(data)])

def compute_radial_distro(data, L=10):
    data = data[-200:, ] # consider final particles position
    
    # # all coordinates are in cartesian coordinates: convert in spherical
    # spherical_data = np.array([np.sqrt(data[i, 0]**2 + data[i, 1]**2 + data[i, 2]**2), # r
    #                            np.arccos(data[i, 2]/np.sqrt(data[i, 0]**2 + data[i, 1]**2 + data[i, 2]**2)), # θ
    #                            np.sign(data[i, 1])*np.arccos(data[i, 0]/np.sqrt(data[i, 0]**2 + data[i, 1]**2)) # φ  
    #                            for i in range(data.shape[0])])
    # count the number of particles at a distance dr
    
    distances = compute_distances(data)
    counter = []
    N, rho = 200, 0.2
    dr = 0.1/2
    eps = 1e-9
    min_dist, max_dist = np.min(distances), np.max(distances)
    # print(f'min_distance = {min_dist}, max_distance = {max_dist}')
    all_r = np.arange(eps, L+eps, dr)
    for r in all_r:
        c = []
        for i, distance_from_p in enumerate(distances):
            # distance_from_p contains the distances of all particles form particle i
            c.append(sum(np.logical_and(distance_from_p > r, distance_from_p < r+dr)))
        counter.append(np.array(c))

    # sum along the rows and normalize 
    radial_distro = np.array(counter).sum(axis=1)/((N-1)*rho*4*np.pi*dr*all_r**2)
    return radial_distro # np.array(counter)#/np.sum(counter)


skip = True
eps = 1e-9
N, rho = 200, 0.2
T_k_simulated = 1.
dt = 10**(-3)
dr = 0.1/2
L = 10.
if not skip:
    for k, sig in enumerate(sigma_cut):
        fig, ax = plt.subplots(1, 3, figsize=(12, 3))
        for j, t in enumerate(types[:3]):
            labels = [r'$V_{LJ}$', r'$T_k^\star$', r'$g(r)$']
            y_labels = [r'$V_{LJ}$', r'$T_k^\star$', r'$g(r)$']
            data = np.loadtxt(f'{folder}{sig}/{t}{ext}')
            time = np.arange(0, len(data), 1)*dt
            if t == 'data':
                data = compute_radial_distro(data)
                print(f'distro integral = {np.trapz(data* 4 * np.pi * rho* np.arange(eps, L+eps, dr)**2, dx=dr)}, for N = {N}, sigma = {sig}')
                print(f'r_min = {np.arange(eps, L+eps, dr)[np.where(data==max(data))]}')
                print(f'Theoretical r_min = {2**(1/6)}')
                mask = np.arange(eps, L+eps, dr) < L/2
                ax[j].plot(np.arange(eps, L+eps, dr)[mask], data[mask], label=labels[j], 
                           color=colors[j], ls='-') # , marker='o' 
                ax[j].set_xlabel(r'$r$')
                ax[j].set_ylim(-0.2, 3.3)
                # ax[j].set_xlim(-0.2, 3.1)
            elif t == 'kinetic_energy':
                # Convert energy into T_k
                T_k = data*2/(3.*N)
                ax[j].plot(time, T_k, label=labels[j], color=colors[j])
                ax[j].hlines(T_k_simulated, ax[j].get_xlim()[0], ax[j].get_xlim()[1], color='k', ls=':')
                ax[j].set_xlabel(r'$t$')
            else:
                ax[j].plot(time, data, label=labels[j], color=colors[j])
                ax[j].set_xlabel(r'$t$')
            ax[j].legend(frameon=False)
            ax[j].grid(ls='--', alpha=0.7)
            ax[j].set_ylabel(y_labels[j])
            # if t != 'data':
            #     ax[j].set_xscale('log')
        plt.tight_layout()
        plt.savefig(f'{saveimg}/energy_{sig}.png', dpi=300)
        plt.close()
        # if sig != sigma_cut[-1]:
        #     plt.show(block=False)
        # else:
        #     plt.show()


plt.rcParams.update({'font.size': 12})  
folder = 'res/Thermostats/'
thermostats = [1, 2] # 
therm_names = ['V_rescaling', 'Andersen']
word = 'kinetic_energy'
T = 2
rho = np.arange(0.08, 0.21, 0.04)   
Vol = 10**3
all_N = np.ceil(rho*Vol).astype(int) 
cols = ['cornflowerblue', 'teal']
skip = True
if not skip:
    for k, therm in enumerate(therm_names):
        print(f'\n{therm_names[k]} thermostat\n')    
        for N in all_N:
            print(f'N = {N}')
            fig, ax = plt.subplots(1, 1, figsize=(12, 3))
            data = np.loadtxt(f'{folder}{therm}/{N}/{word}{ext}')
            data = data[50:]
            time = np.arange(0, len(data), 1)*dt
            ax.plot(time, data, label=r'$E_k$', color=cols[0])
            ax.hlines(data.mean(), ax.get_xlim()[0], ax.get_xlim()[1], label=r'$\langle E_k \rangle$', ls='--', color=cols[1])
            ax.hlines(3/2*N*T, ax.get_xlim()[0], ax.get_xlim()[1], label=r'$\frac{3}{2}NT^\star$', ls=':', color='k')
            ax.set_xlabel(r'$t$')
            ax.set_ylabel(r'$E_k$')
            ax.set_xscale('log')
            ax.legend(frameon=False, ncols=3)
            ax.grid(ls='--', alpha=0.7)
            
            T_k = data*2/(3.*N)
            T_fluctuations = np.var(T_k)/np.mean(T_k)**2

            with open(f'res/Thermostats/{therm}_{N}.dat', 'w') as f:
                f.write(f'Experimental <E_k> = {data.mean()}\n')
                f.write(f'Theoretical <E_k> = {3/2*N*T}\n')

                f.write(f'Experimental T_fluctuations = {T_fluctuations}\n')
                f.write(f'Theoretical T_fluctuations = {2/(3*N)}\n')

            plt.tight_layout()
            plt.savefig(f'{saveimg}/{therm}_{N}.png', dpi=300)
            plt.close()
