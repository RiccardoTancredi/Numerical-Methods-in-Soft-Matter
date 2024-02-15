import numpy as np
import os
import subprocess

def distance(d, t):
    return np.sqrt(np.sum((d-t - np.round((d-t)/10.)*10.)**2, axis=1))

def generate_init_configuration(low=0, high=1, size=(1, 1), sigma=1):
    rows, cols = size
    data = np.random.uniform(low=low, high=high, size=(1, cols))
    for _ in range(1, rows):
        go = True
        while go:
            trial = np.random.uniform(low=low, high=high, size=(1, cols))
            if sum(distance(data, trial) >  sigma) == data.shape[0]:
                data = np.vstack((data, trial))
                go = False

    return data

def read_param():
        
    read_file = np.loadtxt('param.dat', dtype=str)
    var_names, data = [], []
    for ar in read_file:
        for el in ar:
            try:
                data.append(int(float(el)))
            except:
                if '.dat' not in el:
                    var_names.append(el)
                else:
                    if 'velocity' not in el:
                        start_file_name = str(el)
                    else:
                        start_file_name_vel = str(el)
                
    var_names = np.array(var_names)
    var_names = var_names[var_names != '=']
    N, N_strep, box_x, box_y, box_z, T, sigma, sigma_cut = data[:8]

    return N, N_strep, box_x, box_y, box_z, T, sigma, sigma_cut, start_file_name, start_file_name_vel


# Temperatures = ['0.9', '2']  

N, N_strep, box_x, box_y, box_z, T, sigma, sigma_cut, start_file_name, start_file_name_vel = read_param()
file_path_x = f'{start_file_name}'
file_path_v = f'{start_file_name_vel}'

# print(N, N_strep, box_x, box_y, box_z, T, sigma, sigma_cut, start_file_name, start_file_name_vel)

thermostats = [1, 2] # 
therm_names = ['V_rescaling', 'Andersen']
    
T = 2
rho = np.arange(0.08, 0.21, 0.04)   
Vol = box_x*box_y*box_z
all_N = np.ceil(rho*Vol).astype(int) 
for k, therm in enumerate(thermostats):
    print(f'\n{therm_names[k]} thermostat\n')
    for N in all_N:
        print(f'N = {N}')
        os.makedirs(f'res/Thermostats/{therm_names[k]}/{N}', exist_ok=True)
    
        # Update lines in the 'params.dat' file:
        with open('param.dat', 'r') as file:
            lines = file.readlines()
        
        lines[0] = f'N = {N}\n'
        lines[5] = f'temperature = {T}\n'
        lines[12] = f'thermostat = {therm}\n'
        with open('param.dat', 'w') as file:
            file.writelines(lines)

        # Generate initial (x, y, z) positions of the N particles
        init_x = np.random.uniform(low=0, high=box_x, size=(N, 3)) 
        # init_x = generate_init_configuration(low=0, high=box_x, size=(N, 3), sigma=1.2) 
        init_v = np.random.normal(loc=0, scale= np.sqrt(T), size=(N, 3)) 
        np.savetxt(file_path_x, init_x)
        np.savetxt(file_path_v, init_v)

        print('\nFiles created!')

        # Call C program
        executable_file = "mc.out"  

        # Run the compiled executable
        run_command = f"wsl ./{executable_file}"
        run_result = subprocess.run(run_command, shell=True)

        if run_result.returncode != 0:
            print("Error occurred while running the compiled executable.")
            break

with open('param.dat', 'r') as file:
    lines = file.readlines()
# reset initial values
lines[0] = f'N = 200\n' 
lines[5] = f'temperature = 1\n'
lines[12] = f'thermostat = 0\n'
with open('param.dat', 'w') as file:
    file.writelines(lines)
        