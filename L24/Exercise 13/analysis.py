#############
# Cell list #
#############

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

plt.rcParams.update({'font.size': 15})  
folder = 'C-sim/res/'
"""
model = 0: initial algorithm, as the teacher gave us 
           just some code modification: all parameters 
           are provided within the file 'parameters.dat'
model = 1: active dumbbells. Harmonic spring and repulsive 
           force action
model = 2: probe diffusion. Mean square displacement as a 
           function of time. 
"""

model = 2
data_folder = 'OUT'
keywords = ['_active_matter', '_standard', '_f_active_iter_']
files = ['traj', 'x-rel_v']
ext = '.dat'
N = 500 # number of small particels
N_big_particles = 1


colors = ['#780116', '#ff1654', '#c32f27', 
          '#247ba0', '#d8572a', '#db7c26', 
          '#94B0DA', '#70c1b3', '#f7b538', 
          '#b2dbbf']

    
def unwrap_coordinates(coordinates, N, D=2):
    # if coordinates.shape[1] != D:
    #     old_coords = np.copy(coordinates)
    #     coordinates = coordinates[:, :D]

    L = 20
    num_frames = coordinates.shape[0] // N
    # coordinates = coordinates.reshape((num_frames, N, D))
    # unwrapped_coordinates = np.copy(coordinates[:, :D])
    unwrapped_coordinates = coordinates[:, :D] + coordinates[:, D:]*L

    # for i in range(1, num_frames):
    #     displacement = coordinates[i] - coordinates[i-1]
    #     displacement -= np.round(displacement/L)*L
    #     # print(f'i = {i}, displacement = {displacement}\narray: {old_coords[i*N]}')
    #     unwrapped_coordinates[i] = unwrapped_coordinates[i-1] + displacement
    
    return unwrapped_coordinates# .reshape(N*num_frames, D)

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
    msd = np.array([np.mean(np.sum((x_t[i:i+N] - x0)**2, axis=1)) for i in range(0, x_t.shape[0], N)])
    return msd

save_imgs = 'img'
dt = 10**-3
D = 2   # system dimension

all_particles = False
save_fig = False

plt.figure(figsize=(12, 5))
file_to_analyze = 10
for k, num in tqdm(enumerate(range(0, file_to_analyze)), total=file_to_analyze, 
                   desc='Buffering...', colour='green'):
    data_traj = np.loadtxt(f'{folder}/{files[0]}{keywords[2]}{num}{ext}')
    # First row is useless
    # From i = 0 to i = N:
    # Columns: x[i][0], x[i][1], travel[i][0], travel[i][1]
    # Then: x_trap[0], x_trap[1], travel[0][0], travel[0][1]
    # useless row from the next iteration and so on... 
    # (dumbest way of storing data btw)
    offset = 1
    useless_rows = np.array(list(range(0, data_traj.shape[0], 2*offset + (N + N_big_particles))))
    data = np.delete(data_traj, useless_rows, axis=0)
    trap_rows = np.array(list(range(N+N_big_particles, data.shape[0], N+N_big_particles+1)))
    trap_coords = data[trap_rows, :]
    data_all = np.delete(data, trap_rows, axis=0)
    # remove big particle data
    data_big_particle_rows = np.array(list(range(0, data_all.shape[0], N+N_big_particles)))
    data_big_particle = data_all[data_big_particle_rows, :]
    data = np.delete(data_all, data_big_particle_rows, axis=0)

    if all_particles:
        data_unwrapped = unwrap_coordinates(data, N, D)
        msd = mean_squared_displacement(data_unwrapped, N)
    else:
        data_unwrapped = unwrap_coordinates(data_big_particle, 1, D)
        msd = mean_squared_displacement(data_unwrapped, 1)
    
    plt.plot(msd[:3000], label=r'$f_{active}=$'+f'{num}', color=colors[k])

# from scipy.optimize import curve_fit
# popt, pcov = curve_fit(lambda x, a, b: a*x + b, list(range(100)), msd[:100])
# print(popt)
# plt.plot(popt[0]*np.arange(0, 100, 1)+popt[1])
# plt.xlim(0, 200)
plt.xlabel(r'$t$')
plt.ylabel(r'$MSD(t)$')
plt.grid(ls='--', alpha=0.5)
plt.legend(frameon=False, ncols=2) 
plt.tight_layout()
if all_particles and save_fig:
    plt.savefig(f'img/all_particles_diffusion.png', dpi=500)
elif save_fig:
    plt.savefig(f'img/probe_diffusion.png', dpi=500)
plt.show()