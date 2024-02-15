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
    N, N_strep, box_x, box_y, box_z, T, gamma, sigma, mass, sigma_cut = data[:10]

    return N, N_strep, box_x, box_y, box_z, T, gamma, sigma, mass, sigma_cut, start_file_name, start_file_name_vel


# Temperatures = ['0.9', '2']  

N, N_strep, box_x, box_y, box_z, T, gamma, sigma, mass, sigma_cut, start_file_name, start_file_name_vel = read_param()
file_path_x = f'{start_file_name}'
file_path_v = f'{start_file_name_vel}'

# print(N, N_strep, box_x, box_y, box_z, T, sigma, sigma_cut, start_file_name, start_file_name_vel)

Temperatures = np.arange(0.2, 2, 1/5)
Temp_str = [str(t)[:3] for t in Temperatures]

Gamma = list(range(0, 100, 10))
Gamma[0] = 1
Gamma = np.array(Gamma)
Gamma_str = [str(int(g)) for g in Gamma]
    

for integrator in [0, 1]:
    integrator_type = 'Overdamped' if integrator == 0 else 'Underdamped'
    for k, T in enumerate(Temperatures):
        gamma = 1
        os.makedirs(f'res/{integrator_type}/Temperature/{Temp_str[k]}', exist_ok=True)
        
        # Update lines in the 'params.dat' file:
        with open('param.dat', 'r') as file:
            lines = file.readlines()
        lines[15] = f'integrator = {integrator}\n'
        lines[5] = f'temperature = {T}\n'
        lines[6] = f'gamma = {gamma}\n'
        lines[14] = f'model = {0}\n'
        with open('param.dat', 'w') as file:
            file.writelines(lines)

        # Generate initial (x, y, z) positions of the N particles
        # init_x = np.random.uniform(low=0, high=box_x, size=(N, 3)) 
        init_x = generate_init_configuration(low=0, high=box_x, size=(N, 3), sigma=1.2) 
        init_v = np.random.normal(loc=0, scale= 1., size=(N, 3)) 
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


    for l, gamma in enumerate(Gamma):
        T = 1
        os.makedirs(f'res/{integrator_type}/Gamma/{Gamma_str[l]}', exist_ok=True)
        
        # Update lines in the 'params.dat' file:
        with open('param.dat', 'r') as file:
            lines = file.readlines()
        lines[15] = f'integrator = {integrator}\n'
        lines[5] = f'temperature = {T}\n'
        lines[6] = f'gamma = {gamma}\n'
        lines[14] = f'model = {1}\n'
        with open('param.dat', 'w') as file:
            file.writelines(lines)
        
        # Generate initial (x, y, z) positions of the N particles
        # init_x = np.random.uniform(low=0, high=box_x, size=(N, 3)) 
        init_x = generate_init_configuration(low=0, high=box_x, size=(N, 3), sigma=1.2) 
        init_v = np.random.normal(loc=0, scale= 1., size=(N, 3)) 
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
lines[5] = f'temperature = {1}\n'
lines[6] = f'gamma = {1}\n'
lines[14] = f'model = {0}\n'
lines[15] = f'integrator = {0}\n' # reset 
with open('param.dat', 'w') as file:
    file.writelines(lines)
        