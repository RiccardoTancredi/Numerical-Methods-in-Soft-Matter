import numpy as np
import matplotlib.pyplot as plt
import os
from tqdm import tqdm

plt.rcParams.update({'font.size': 14})  
folder = 'SOURCE/res/Ex2'

N = 200

Temperatures = np.arange(0.2, 2, 1/5)
Temp_str = [str(t)[:3] for t in Temperatures]

Gamma = list(range(10, 100, 10))
Gamma = np.array(Gamma)
Gamma_str = [str(int(g)) for g in Gamma]
    
Kappa = list(range(1, 10, 1))
Kappa_str = [str(K) for K in Kappa]

integrator_type = 'Overdamped'

spring = [[0, 0], [2, 0]]

file_type = ['data', 'energy', 'kinetic_energy', 'velocity']
ext = '.dat'

colors = ["#247ba0", "#70c1b3", "#b2dbbf",
          "#ff1654", "#780116", "#f7b538", 
          "#db7c26", "#d8572a", "#c32f27"]

    
def unwrap_coordinates(coordinates, N):
    L = 20
    num_frames = coordinates.shape[0] // N
    coordinates = coordinates.reshape((num_frames, N, 3))
    unwrapped_coordinates = np.copy(coordinates)
    num_frames, num_particles, _ = coordinates.shape
    
    for i in range(1, num_frames):
        displacement = coordinates[i] - coordinates[i-1]
        displacement -= np.round(displacement/L)*L
        unwrapped_coordinates[i] = unwrapped_coordinates[i-1] + displacement
    
    return unwrapped_coordinates.reshape(N*num_frames, 3)

def pbc(x1, x2):
    L = 20
    return x1-x2 - np.round((x1-x2)/L)*L
    
    
def mean_squared_displacement(matrix, N):
    """
    Args:
        matrix (np.array): data matrix containing the positions
        N (int): number of small particles
    """
    x0 = matrix[:N, :]
    x_t = matrix[N:, :]
    msd= np.array([np.mean(np.sum((x_t[i:i+N]- x0)**2, axis=1)) for i in range(0, x_t.shape[0], N)])
    return msd

save_imgs = 'SOURCE/img/Ex2'
dt = 10**-3

def correlation(velocity, N):
    v0 = velocity[0:N, :]
    vt = velocity[N:, :]
    t_max = (velocity.shape[0])//N
    correlation_function = np.zeros(shape=t_max-1)
    
    for t in range(t_max-1):
        correlation_function[t] = np.mean(np.sum(velocity[t*N:(t+1)*N] * v0, axis=1))
        # correlation_function[t] = np.mean(np.sum(velocity[t*N:] * velocity[:(t_max-t)*N], axis=1))
    return correlation_function #  / correlation_function[0]

def correlation_time(data, N, dt):
    t_max = data.shape[0]
    all_tau_x, all_tau_y = [], []
    for n in range(N):
        correlation_function_x = np.zeros(shape=t_max)
        correlation_function_y = np.zeros(shape=t_max)
        for t in range(t_max):
            correlation_function_x[t] = np.sum(data[:t_max-t, n, 0]*data[t:t_max, n, 0])/(t_max-t)-((1./(t_max-t))**2)*np.sum(data[0:t_max-t, n, 0])*np.sum(data[t:t_max, n, 0])
            correlation_function_y[t] = np.sum(data[:t_max-t, n, 1]*data[t:t_max, n, 1])/(t_max-t)-((1./(t_max-t))**2)*np.sum(data[0:t_max-t, n, 1])*np.sum(data[t:t_max, n, 1])
        correlation_function_x /= correlation_function_x[0]
        correlation_function_y /= correlation_function_y[0]
        
        # Now that we have our correlation function: -> estimate t_O = integrated autocorrelation time
        wh_x = np.where(correlation_function_x < 0)[0][0]
        wh_y = np.where(correlation_function_y < 0)[0][0]
        tau_x = np.trapz(correlation_function_x[:wh_x], dx=dt)
        tau_y = np.trapz(correlation_function_y[:wh_y], dx=dt)
        # print('\n', tau_x, tau_y)
        # plt.plot(correlation_function_x[:wh_x], label='x')
        # plt.plot(correlation_function_y[:wh_y], label='y')
        # plt.legend()
        # plt.show()
        all_tau_x.append(tau_x)
        all_tau_y.append(tau_y)
    
    return np.array(all_tau_x), np.array(all_tau_y)

skip = True

dim = 2     # system's dimension
c_t_folder = 'SOURCE/corr_time'

if not skip:
    
    ###############
    # Temperature #
    ###############

    all_means_x, all_vars_x = [], []
    all_means_y, all_vars_y = [], []
    corr_time_x, corr_time_y = [], []
    corr_time_var_x, corr_time_var_y = [], []
    for k, T in tqdm(enumerate(Temperatures), total=len(Temperatures), desc='Processing', colour='green'):
        gamma = 1
        K = 1
        data_path = f'{folder}/{integrator_type}/Temperature/{Temp_str[k]}/{file_type[0]}{ext}'
        data = np.loadtxt(data_path)
        data = unwrap_coordinates(data, N)[:, :dim] # z coordinate has been used
        # data is an array containing the positions along 
        # x, and y of the N = 200 simulated particles
        # The idea is to compute the average for each particle individually 
        # and then average among the N particles 
        num_frames = data.shape[0] // N
        data_reshaped = data.reshape((num_frames, N, dim))
        mean_x = np.mean([data_reshaped[:, particle, 0].mean() for particle in range(N)])
        var_x = np.var([data_reshaped[:, particle, 0].mean() for particle in range(N)])/N
        mean_y = np.mean([data_reshaped[:, particle, 1].mean() for particle in range(N)])
        var_y = np.var([data_reshaped[:, particle, 1].mean() for particle in range(N)])/N
        all_means_x.append(mean_x)
        all_vars_x.append(var_x)
        all_means_y.append(mean_y)
        all_vars_y.append(var_y)
        taus = correlation_time(data_reshaped, N, dt)
        corr_time_x.append(np.mean(taus[0]))
        corr_time_y.append(np.mean(taus[1]))
        corr_time_var_x.append(np.var(taus[0])/N)
        corr_time_var_y.append(np.var(taus[1])/N)
    np.savetxt(f'{c_t_folder}/Temp_mean_x.dat', corr_time_x)
    np.savetxt(f'{c_t_folder}/Temp_vars_x.dat', corr_time_var_x)
    np.savetxt(f'{c_t_folder}/Temp_mean_y.dat', corr_time_y)
    np.savetxt(f'{c_t_folder}/Temp_vars_y.dat', corr_time_var_y)

    
    # Means and Stds
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(15, 5)) 
    ax[0].errorbar(Temperatures, all_means_x, np.sqrt(all_vars_x), color=colors[0], fmt='o',
                   label=r'$K\sigma^2/\varepsilon\equiv 1$'+f'\n'+r'$\gamma\tau\equiv 1$')
    ax[0].set_xlabel(r'$T^\star$')
    ax[0].set_ylabel(r'$\overline{x}$')
    ax[0].set_xticks(Temperatures, Temp_str)
    ax[0].legend(frameon=False, loc='upper left')
    ax[0].grid(ls='--', alpha=0.5)
    
    ax[1].errorbar(Temperatures, all_means_y, np.sqrt(all_vars_y), color=colors[-1], fmt='o',
                   label=r'$K\sigma^2/\varepsilon\equiv 1$'+f'\n'+r'$\gamma\tau\equiv 1$')
    ax[1].set_xlabel(r'$T^\star$')
    ax[1].set_ylabel(r'$\overline{y}$')
    ax[1].set_xticks(Temperatures, Temp_str)
    ax[1].legend(frameon=False, loc='upper left')
    ax[1].grid(ls='--', alpha=0.5)

    plt.tight_layout()
    plt.savefig(f'{save_imgs}/Temperature.png', dpi=500)
    plt.close()

    
    # Correlation time
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(15, 5)) 
    ax[0].errorbar(Temperatures, corr_time_x, np.sqrt(corr_time_var_x), color=colors[1], fmt='o',
                   label=r'$K\sigma^2/\varepsilon\equiv 1$'+f'\n'+r'$\gamma\tau\equiv 1$')
    ax[0].set_xlabel(r'$T^\star$')
    ax[0].set_ylabel(r'$\overline{\tau}_x$')
    ax[0].set_xticks(Temperatures, Temp_str)
    ax[0].legend(frameon=False, loc='upper left')
    ax[0].grid(ls='--', alpha=0.5)
    
    ax[1].errorbar(Temperatures, corr_time_y, np.sqrt(corr_time_var_y), color=colors[-2], fmt='o',
                   label=r'$K\sigma^2/\varepsilon\equiv 1$'+f'\n'+r'$\gamma\tau\equiv 1$')
    ax[1].set_xlabel(r'$T^\star$')
    ax[1].set_ylabel(r'$\overline{\tau}_y$')
    ax[1].set_xticks(Temperatures, Temp_str)
    ax[1].legend(frameon=False, loc='upper left')
    ax[1].grid(ls='--', alpha=0.5)

    plt.tight_layout()
    plt.savefig(f'{save_imgs}/Temperature_Corr_Time.png', dpi=500)
    plt.close()


    #########
    # Gamma #
    #########

    all_means_x, all_vars_x = [], []
    all_means_y, all_vars_y = [], []
    corr_time_x, corr_time_y = [], []
    corr_time_var_x, corr_time_var_y = [], []
    for l, gamma in tqdm(enumerate(Gamma), total=len(Gamma), desc='Processing', colour='green'):
        T = 1
        K = 1
        data_path = f'{folder}/{integrator_type}/Gamma/{Gamma_str[l]}/{file_type[0]}{ext}'
        data = np.loadtxt(data_path)
        data = unwrap_coordinates(data, N)[:, :dim] # z coordinate has been used
        num_frames = data.shape[0] // N
        data_reshaped = data.reshape((num_frames, N, dim))
        mean_x = np.mean([data_reshaped[:, particle, 0].mean() for particle in range(N)])
        var_x = np.var([data_reshaped[:, particle, 0].mean() for particle in range(N)])/N
        mean_y = np.mean([data_reshaped[:, particle, 1].mean() for particle in range(N)])
        var_y = np.var([data_reshaped[:, particle, 1].mean() for particle in range(N)])/N
        all_means_x.append(mean_x)
        all_vars_x.append(var_x)
        all_means_y.append(mean_y)
        all_vars_y.append(var_y)
        taus = correlation_time(data_reshaped, N, dt)
        corr_time_x.append(np.mean(taus[0]))
        corr_time_y.append(np.mean(taus[1]))
        corr_time_var_x.append(np.var(taus[0])/N)
        corr_time_var_y.append(np.var(taus[1])/N)
    np.savetxt(f'{c_t_folder}/Gamma_mean_x.dat', corr_time_x)
    np.savetxt(f'{c_t_folder}/Gamma_vars_x.dat', corr_time_var_x)
    np.savetxt(f'{c_t_folder}/Gamma_mean_y.dat', corr_time_y)
    np.savetxt(f'{c_t_folder}/Gamma_vars_y.dat', corr_time_var_y)

    
    # Means and Stds
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(15, 5)) 
    ax[0].errorbar(Gamma, all_means_x, np.sqrt(all_vars_x), color=colors[0], fmt='o',
                   label=r'$T^\star\equiv 1$'+f'\n'+r'$K\sigma^2/\varepsilon\equiv 1$')
    ax[0].set_xlabel(r'$\gamma\tau$')
    ax[0].set_ylabel(r'$\overline{x}$')
    ax[0].set_xticks(Gamma, Gamma_str)
    ax[0].legend(frameon=False)
    ax[0].grid(ls='--', alpha=0.5)
    
    ax[1].errorbar(Gamma, all_means_y, np.sqrt(all_vars_y), color=colors[-1], fmt='o',
                   label=r'$T^\star\equiv 1$'+f'\n'+r'$K\sigma^2/\varepsilon\equiv 1$')
    ax[1].set_xlabel(r'$\gamma\tau$')
    ax[1].set_ylabel(r'$\overline{y}$')
    ax[1].set_xticks(Gamma, Gamma_str)
    ax[1].legend(frameon=False)
    ax[1].grid(ls='--', alpha=0.5)

    plt.tight_layout()
    plt.savefig(f'{save_imgs}/Gamma.png', dpi=500)
    plt.close()
    
    # Correlation time
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(15, 5)) 
    ax[0].errorbar(Gamma, corr_time_x, np.sqrt(corr_time_var_x), color=colors[1], fmt='o',
                   label=r'$T^\star\equiv 1$'+f'\n'+r'$K\sigma^2/\varepsilon\equiv 1$')
    ax[0].set_xlabel(r'$\gamma\tau$')
    ax[0].set_ylabel(r'$\overline{\tau}_x$')
    ax[0].set_xticks(Gamma, Gamma_str)
    ax[0].legend(frameon=False)
    ax[0].grid(ls='--', alpha=0.5)
    
    ax[1].errorbar(Gamma, corr_time_y, np.sqrt(corr_time_var_y), color=colors[-2], fmt='o',
                   label=r'$T^\star\equiv 1$'+f'\n'+r'$K\sigma^2/\varepsilon\equiv 1$')
    ax[1].set_xlabel(r'$\gamma\tau$')
    ax[1].set_ylabel(r'$\overline{\tau}_y$')
    ax[1].set_xticks(Gamma, Gamma_str)
    ax[1].legend(frameon=False)
    ax[1].grid(ls='--', alpha=0.5)

    plt.tight_layout()
    plt.savefig(f'{save_imgs}/Gamma_Corr_Time.png', dpi=500)
    plt.close()
    
    #########
    # Kappa #
    #########

    all_means_x, all_vars_x = [], []
    all_means_y, all_vars_y = [], []
    corr_time_x, corr_time_y = [], []
    corr_time_x, corr_time_y = [], []
    corr_time_var_x, corr_time_var_y = [], []
    for j, K in tqdm(enumerate(Kappa), total=len(Kappa), desc='Processing', colour='green'):
        T = 1
        gamma = 1
        data_path = f'{folder}/{integrator_type}/Kappa/{Kappa_str[j]}/{file_type[0]}{ext}'
        data = np.loadtxt(data_path)
        data = unwrap_coordinates(data, N)[:, :dim] # z coordinate has been used
        num_frames = data.shape[0] // N
        data_reshaped = data.reshape((num_frames, N, dim))
        mean_x = np.mean([data_reshaped[:, particle, 0].mean() for particle in range(N)])
        var_x = np.var([data_reshaped[:, particle, 0].mean() for particle in range(N)])/N
        mean_y = np.mean([data_reshaped[:, particle, 1].mean() for particle in range(N)])
        var_y = np.var([data_reshaped[:, particle, 1].mean() for particle in range(N)])/N
        all_means_x.append(mean_x)
        all_vars_x.append(var_x)
        all_means_y.append(mean_y)
        all_vars_y.append(var_y)
        taus = correlation_time(data_reshaped, N, dt)
        corr_time_x.append(np.mean(taus[0]))
        corr_time_y.append(np.mean(taus[1]))
        corr_time_var_x.append(np.var(taus[0])/N)
        corr_time_var_y.append(np.var(taus[1])/N)
    np.savetxt(f'{c_t_folder}/Kappa_mean_x.dat', corr_time_x)
    np.savetxt(f'{c_t_folder}/Kappa_vars_x.dat', corr_time_var_x)
    np.savetxt(f'{c_t_folder}/Kappa_mean_y.dat', corr_time_y)
    np.savetxt(f'{c_t_folder}/Kappa_vars_y.dat', corr_time_var_y)

        
    # Means and Stds    
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(15, 5)) 
    ax[0].errorbar(Kappa, all_means_x, np.sqrt(all_vars_x), color=colors[0], fmt='o',
                   label=r'$T^\star\equiv 1$'+f'\n'+r'$\gamma\tau\equiv 1$')
    ax[0].set_xlabel(r'$K\sigma^2/\varepsilon$')
    ax[0].set_ylabel(r'$\overline{x}$')
    ax[0].set_xticks(Kappa, Kappa_str)
    ax[0].legend(frameon=False)
    ax[0].grid(ls='--', alpha=0.5)
    
    ax[1].errorbar(Kappa, all_means_y, np.sqrt(all_vars_y), color=colors[-1], fmt='o',
                   label=r'$T^\star\equiv 1$'+f'\n'+r'$\gamma\tau\equiv 1$')
    ax[1].set_xlabel(r'$K\sigma^2/\varepsilon$')
    ax[1].set_ylabel(r'$\overline{y}$')
    ax[1].set_xticks(Kappa, Kappa_str)
    ax[1].legend(frameon=False)
    ax[1].grid(ls='--', alpha=0.5)

    plt.tight_layout()
    plt.savefig(f'{save_imgs}/Kappa.png', dpi=500)
    plt.close()
    
    # Correlation time    
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(15, 5)) 
    ax[0].errorbar(Kappa, corr_time_x, np.sqrt(corr_time_var_y), color=colors[1], fmt='o',
                   label=r'$T^\star\equiv 1$'+f'\n'+r'$\gamma\tau\equiv 1$')
    ax[0].set_xlabel(r'$K\sigma^2/\varepsilon$')
    ax[0].set_ylabel(r'$\overline{\tau}_x$')
    ax[0].set_xticks(Kappa, Kappa_str)
    ax[0].legend(frameon=False)
    ax[0].grid(ls='--', alpha=0.5)
    
    ax[1].errorbar(Kappa, corr_time_y, np.sqrt(corr_time_var_y), color=colors[-2], fmt='o',
                   label=r'$T^\star\equiv 1$'+f'\n'+r'$\gamma\tau\equiv 1$')
    ax[1].set_xlabel(r'$K\sigma^2/\varepsilon$')
    ax[1].set_ylabel(r'$\overline{\tau}_y$')
    ax[1].set_xticks(Kappa, Kappa_str)
    ax[1].legend(frameon=False)
    ax[1].grid(ls='--', alpha=0.5)

    plt.tight_layout()
    plt.savefig(f'{save_imgs}/Kappa_Corr_Time.png', dpi=500)
    plt.close()

#############
# 2 springs #
#############

skip = False
initials = [[0, 0], [0.125, 0.125], [0.25, 0.25], [0.5, 0.5], [0.75, 0.75], [1, 1]]
initials_str = ['0.0',  '0.1', '0.2', '0.5', '0.8', '1.0']
initials_plot = [r'$(0, 0)^T$', r'$(0.125, 0.125)^T$', r'$(0.25, 0.25)^T$',
                 r'$(0.5, 0.5)^T$', r'$(0.75, 0.75)^T$', r'$(1, 1)^T$']
initials_plot = [r'$0$', r'$0.125$', r'$0.25$', r'$0.5$', r'$0.75$', r'$1$']

spring = np.array([[0, 0], [2, 0]])
if not skip:
    max_if_not_found = 0 # 1./dt * 10**2
    first_passage_time_to_1_mean = []
    first_passage_time_to_1_var = []
    first_passage_time_to_2_mean = []
    first_passage_time_to_2_var = []
    for index, initial in tqdm(enumerate(initials), total=len(initials_str), desc='Loading...', colour='green'):
        data_path = f'{folder}/{integrator_type}/2_springs/{initials_str[index]}/{file_type[0]}{ext}'
        data = np.loadtxt(data_path)
        # data = unwrap_coordinates(data, N)[:, :dim] # z coordinate has been used
        data = data[:, :dim]
        num_frames = data.shape[0] // N
        # data_reshaped = data.reshape((num_frames, N, dim))
        distance_from_spring_1 = pbc(data, spring[0])
        distance_from_spring_2 = pbc(data, spring[1])
        distance_1 = np.sqrt(np.sum(distance_from_spring_1**2, axis=1)).reshape(num_frames, N)
        distance_2 = np.sqrt(np.sum(distance_from_spring_2**2, axis=1)).reshape(num_frames, N)

        first_passage_time_to_1, first_passage_time_to_2 = [], []
        for n in range(N):
            try:
                # FPT_1 = np.where((distance_1[:, n] < distance_2[:, n]) == True)[0][0]
                FPT_2 = np.where((distance_2[:, n] < distance_1[:, n]) == True)[0][0]
            except:
                # FPT_1 = max_if_not_found
                FPT_2 = max_if_not_found
                print('Unable to find FPT')
            # first_passage_time_to_1.append(FPT_1)
            first_passage_time_to_2.append(FPT_2)
        
        # first_passage_time_to_1_mean.append(np.mean(first_passage_time_to_1))
        # first_passage_time_to_1_var.append(np.var(first_passage_time_to_1))
        first_passage_time_to_2_mean.append(np.mean(first_passage_time_to_2))
        first_passage_time_to_2_var.append(np.var(first_passage_time_to_2)/len(first_passage_time_to_2))

    # print(f'first_passage_time_to_1_mean = {first_passage_time_to_1_mean}')
    # print(f'first_passage_time_to_1_var = {np.sqrt(first_passage_time_to_1_var/N)}')
    print(f'first_passage_time_to_2_mean = {first_passage_time_to_2_mean}')
    print(f'first_passage_time_to_2_var = {first_passage_time_to_2_var}')

    # Plot of first mean passage time to jump in the trap of the second spring
    # as a function of the initial position
    plt.errorbar([initials[i][0] for i in range(len(initials_str))], 
                 np.array(first_passage_time_to_2_mean)*dt, 
                 np.sqrt(first_passage_time_to_2_var)*dt, fmt='o', ls='--',
                 color=colors[3], label=r'$FPT \quad 1\to 2$')
    plt.xlabel(r'$x$') # \vec{r}(t_0) 
    plt.xticks([initials[i][0] for i in range(len(initials_str))], initials_plot)
    plt.ylabel(r'$FPT\quad 1\to 2$') 
    # plt.yticks(np.array(first_passage_time_to_2_mean), np.array(first_passage_time_to_2_mean))
    # plt.yscale('log')
    plt.legend()
    plt.grid(ls='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(f'{save_imgs}/2_springs_FPT.png', dpi=500)
    plt.show()