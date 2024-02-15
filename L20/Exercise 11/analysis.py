import numpy as np
import matplotlib.pyplot as plt
import os
from tqdm import tqdm

plt.rcParams.update({'font.size': 14})  
folder = 'SOURCE/res'

N = 200

Temperatures = np.arange(0.2, 2, 1/5)
Temp_str = [str(t)[:3] for t in Temperatures]

Gamma = list(range(10, 100, 10))
Gamma = np.array(Gamma)
Gamma_str = [str(int(g)) for g in Gamma]
    
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

save_imgs = 'SOURCE/img'
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

skip = False

if not skip:
    for integrator in tqdm([0, 1]):
        diff_coeffs = []
        integrator_type = 'Overdamped' if integrator == 0 else 'Underdamped'
        plt.figure(figsize=(12, 5))
        for k, T in enumerate(Temperatures):
            gamma = 1
            data_path = f'{folder}/{integrator_type}/Temperature/{Temp_str[k]}/{file_type[0]}{ext}'
            data_wrapped = np.loadtxt(data_path)
            data = unwrap_coordinates(data_wrapped, N)
            MSD = mean_squared_displacement(data, N)
            time = np.linspace(0, len(MSD)*dt, len(MSD))
            
            if integrator == 1:
                # Compute diffusion coefficient
                velocity = np.loadtxt(f'{folder}/{integrator_type}/Temperature/{Temp_str[k]}/{file_type[-1]}{ext}')
                en_k = np.loadtxt(f'{folder}/{integrator_type}/Temperature/{Temp_str[k]}/{file_type[2]}{ext}')
                wh = np.where(en_k < en_k.mean())[0][0] if en_k[0] > en_k.mean() else 0
                correlation_func = correlation(velocity, N)
                # np.savetxt(f'SOURCE/VACF/veloc_T_{k}.dat', correlation_func)
                # wh = np.where(correlation_func < 0)[0][0]
                D = np.trapz(correlation_func[wh:], dx=dt)/3.
                diff_coeffs.append(D)

            plt.plot(time, MSD, label=r'$T^\star=$'+f'{Temp_str[k]}', color=colors[k])
        plt.xlabel(r'$t$')
        plt.ylabel(r'$MSD(t)$')
        plt.legend(frameon=False, ncols=3)
        plt.grid(ls='--', alpha=0.5)
        plt.tight_layout()
        plt.savefig(f'{save_imgs}/{integrator_type}/Temperature.png', dpi=500)
        plt.close()

        if integrator == 1:
            np.savetxt(f'SOURCE/Diff_Temp.dat', np.array(diff_coeffs))
            
        diff_coeffs = []
        plt.figure(figsize=(12, 5))
        for l, gamma in enumerate(Gamma):
            T = 1
            data_path = f'{folder}/{integrator_type}/Gamma/{Gamma_str[l]}/{file_type[0]}{ext}'
            data_wrapped = np.loadtxt(data_path)
            data = unwrap_coordinates(data_wrapped, N)
            MSD = mean_squared_displacement(data, N)
            time = np.linspace(0, len(MSD)*dt, len(MSD))
            
            if integrator == 1:
                # Compute diffusion coefficient
                en_k = np.loadtxt(f'{folder}/{integrator_type}/Gamma/{Gamma_str[l]}/{file_type[2]}{ext}')
                # wh = np.where(en_k < en_k.mean())[0][0] if en_k[0] > en_k.mean() else 0
                velocity = np.loadtxt(f'{folder}/{integrator_type}/Gamma/{Gamma_str[l]}/{file_type[-1]}{ext}')
                correlation_func = correlation(velocity, N)
                # np.savetxt(f'SOURCE/VACF/veloc_gamma_{l}.dat', correlation_func)
                wh = np.where(correlation_func < 0)[0][0]
                D = np.trapz(correlation_func[:wh], dx=dt)/3.
                diff_coeffs.append(D)

            
            plt.plot(time, MSD, label=r'$\gamma\tau=$'+f'{Gamma_str[l]}', color=colors[l])
        plt.xlabel(r'$t$')
        plt.ylabel(r'$MSD(t)$')
        plt.legend(frameon=False, ncols=3)
        plt.grid(ls='--', alpha=0.5)
        plt.tight_layout()
        plt.savefig(f'{save_imgs}/{integrator_type}/Gamma.png', dpi=500)
        plt.close()

        if integrator == 1:
            np.savetxt(f'SOURCE/Diff_Gamma.dat', np.array(diff_coeffs))
            

skip = False

if not skip:
    diff_T = np.loadtxt('SOURCE/Diff_Temp.dat')
    plt.errorbar(Temperatures, diff_T, np.std(diff_T), fmt='o', 
                 color='royalblue', label=r'$\gamma\tau \equiv 1$')
    # plt.scatter(Temperatures, diff_T, color='royalblue', 
    #             label=r'$\gamma\tau \equiv 1$', s=50)
    # lims = plt.xlim()
    # plt.hlines(np.mean(diff_T), *lims, color='k', ls='--', label=r'$\overline{D}(T^\star) = $'+f'{round(np.mean(diff_T), 3)}')
    # plt.xlim(lims)
    plt.xlabel(r'$T^\star$')
    plt.xticks(Temperatures, Temp_str)
    plt.ylabel(r'$D(T^\star)$')
    plt.legend(loc='lower right', frameon=False, ncols=1)
    plt.grid(ls='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(f'{save_imgs}/Diff_Temp.png', dpi=500)
    plt.close()

    diff_gamma = np.loadtxt('SOURCE/Diff_Gamma.dat')
    plt.errorbar(Gamma, diff_gamma, np.std(diff_gamma), fmt='o', 
                 color='cornflowerblue', label=r'$T^\star \equiv 1$')
    # plt.scatter(Gamma, diff_gamma, color='cornflowerblue', 
    #             label=r'$T^\star \equiv 1$', s=50)
    # lims = plt.xlim()
    # plt.hlines(np.mean(diff_gamma), *lims, color='k', ls='--', 
    #            label=r'$\overline{D}(\gamma\tau) = $'+f'{round(np.mean(diff_gamma), 3)}')
    # plt.xlim(lims)
    plt.xlabel(r'$\gamma\tau$')
    plt.xticks(Gamma, Gamma_str)
    plt.ylabel(r'$D(\gamma\tau)$')
    plt.legend(frameon=False, ncols=1)
    plt.grid(ls='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(f'{save_imgs}/Diff_Gamma.png', dpi=500)
    plt.close()

