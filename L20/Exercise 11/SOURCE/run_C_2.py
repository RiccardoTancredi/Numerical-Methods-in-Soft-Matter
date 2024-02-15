import numpy as np
import os
import subprocess

def distance(d, t):
    return np.sqrt(np.sum((d-t - np.round((d-t)/10.)*10.)**2, axis=1))

def generate_init_configuration(low=0, high=1, size=(1, 1), n_spring=1):
    rows, cols = size
    # all particles starts from the same initial position
    if n_spring == 2:
        data = np.array([low, high, 0]*rows)
    else:
        data = np.array([1, 1, 0] * rows)
    data = data.reshape((rows, cols))
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
    
Kappa = list(range(1, 10, 1))
Kappa_str = [str(K) for K in Kappa]

integrator_type = 'Overdamped'

num_spring = 0
skip = True
if not skip:
    for k, T in enumerate(Temperatures):
        gamma = 1
        K = 1
        os.makedirs(f'res/Ex2/{integrator_type}/Temperature/{Temp_str[k]}', exist_ok=True)
        
        # Update lines in the 'params.dat' file:
        with open('param.dat', 'r') as file:
            lines = file.readlines()
        lines[-1] = f'exercise = 2\n'
        lines[7] = f'K = {K}\n'
        lines[5] = f'temperature = {T}\n'
        lines[6] = f'gamma = {gamma}\n'
        lines[15] = f'model = {0}\n'
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
        K = 1
        os.makedirs(f'res/Ex2/{integrator_type}/Gamma/{Gamma_str[l]}', exist_ok=True)
        
        # Update lines in the 'params.dat' file:
        with open('param.dat', 'r') as file:
            lines = file.readlines()
        lines[-1] = f'exercise = 2\n'
        lines[7] = f'K = {K}\n'
        lines[5] = f'temperature = {T}\n'
        lines[6] = f'gamma = {gamma}\n'
        lines[15] = f'model = {1}\n'
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

    for j, K in enumerate(Kappa):
        T = 1
        gamma = 1
        os.makedirs(f'res/Ex2/{integrator_type}/Kappa/{Kappa_str[j]}', exist_ok=True)
        
        # Update lines in the 'params.dat' file:
        with open('param.dat', 'r') as file:
            lines = file.readlines()
        lines[-1] = f'exercise = 2\n'
        lines[7] = f'K = {K}\n'
        lines[5] = f'temperature = {T}\n'
        lines[6] = f'gamma = {gamma}\n'
        lines[15] = f'model = {2}\n'
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

skip = False
initials = [[0, 0], [0.125, 0.125], [0.25, 0.25], [0.5, 0.5], [0.75, 0.75], [1, 1]]
initials_str = ['0.0',  '0.1', '0.2', '0.5', '0.8', '1.0']

spring = [[0, 0], [2, 0]]
num_spring = 1
if not skip:
    gamma, T, K = [1, 1, 1]
    for h, initial in enumerate(initials):
        os.makedirs(f'res/Ex2/{integrator_type}/2_springs/{initials_str[h]}', exist_ok=True)
            
        # Update lines in the 'params.dat' file:
        with open('param.dat', 'r') as file:
            lines = file.readlines()
        lines[-1] = f'exercise = 2\n'
        lines[7] = f'K = {K}\n'
        lines[5] = f'temperature = {T}\n'
        lines[6] = f'gamma = {gamma}\n'
        lines[15] = f'model = {3}\n'
        lines[19] = f'spring_2_x = {spring[num_spring][0]}\n'
        lines[20] = f'spring_2_y = {spring[num_spring][1]}\n'
        with open('param.dat', 'w') as file:
            file.writelines(lines)

        init_x = generate_init_configuration(low=initial[0], high=initial[1], size=(N, 3), n_spring=2) 
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

# reset
with open('param.dat', 'r') as file:
    lines = file.readlines()
lines[-1] = f'exercise = 1\n'
lines[7] = f'K = {1}\n'
lines[5] = f'temperature = {1}\n'
lines[6] = f'gamma = {1}\n'
lines[15] = f'model = {0}\n'
lines[19] = f'spring_2_x = {spring[num_spring-1][0]}\n'
lines[20] = f'spring_2_y = {spring[num_spring-1][1]}\n'
    
with open('param.dat', 'w') as file:
    file.writelines(lines)
        