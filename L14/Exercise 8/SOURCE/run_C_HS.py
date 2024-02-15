import numpy as np
import os
import subprocess

def generate_init_configuration(low=0, high=1, size=(1, 1), sigma=1):
    # accelerate converge by putting all points at a distance sigma
    rows, cols = size
    data = np.random.uniform(low=low, high=high, size=(1, cols))
    for _ in range(1, rows):
        go = True
        while go:
            trial = np.random.uniform(low=low, high=high, size=(1, cols))
            if sum(np.sqrt(((data-trial)**2).sum(axis=1)) >  sigma) == data.shape[0]:
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
                    
    var_names = np.array(var_names)
    var_names = var_names[var_names != '=']
    keyword = var_names[np.where(var_names == 'kind')[0][0]+1]             
                    
    N, N_strep, box_x, box_y, box_z, T, sigma, delta_x = data[:8]
    iteration = data[-1]
    with open('param.dat', 'r') as file:
        lines = file.readlines()
        delta_x = float(lines[7].split('=')[-1])
    delta_x = delta_x if delta_x != 1. else int(1) 

    return N, N_strep, box_x, box_y, box_z, T, sigma, delta_x, iteration, str(keyword), start_file_name

N, N_strep, box_x, box_y, box_z, T, sigma, delta_x, iteration, keyword, start_file_name = read_param()

Volumes = np.array([5, 6, 7, 13, 22])
all_delta_x = ['0.01', '0.05', '0.3', '0.5', '1']
file_path_x = f'{start_file_name}'

for ty in ["random", "cube"][::-1]:
    run = 0
    while run < len(Volumes):
        V = Volumes[run]
        print(f'\nBuffering V = {V}\n')
        
        for d in all_delta_x:
            with open('param.dat', 'r') as file:
                lines = file.readlines()
            
            my_dict = {2:'x', 3:'y', 4:'z'}
            for i in range(2, 5):
                lines[i] = f'box_{my_dict[i]} = {V}\n' 
            lines[7] = f'delta_x = {d}\n'    # update d_max

            lines[-2] = f'myseed = 02012024\n'  # reset seed    
            lines[-1] = f'iteration = 1\n'      # reset counter
            with open('param.dat', 'w') as file:
                file.writelines(lines)

            it = 1
            while it < 11:
                print(f'\nBuffering iteration = {it}')
                
                N, N_strep, box_x, box_y, box_z, T, sigma, delta_x, iteration, keyword, start_file_name = read_param()
                # Generate initial (x, y, z) positions of the N particles
                init_x = np.random.uniform(low=0, high=box_x, size=(N, 3))
                np.savetxt(file_path_x, init_x)
                print('\nFile created!')

                # Create folder to store C output results
                os.makedirs(f'res/{keyword}/{box_x}_{delta_x}/{iteration}', exist_ok=True)

                # Call C program
                executable_file = "mc.out"  

                # Run the compiled executable
                run_command = f"wsl ./{executable_file}"
                run_result = subprocess.run(run_command, shell=True)

                if run_result.returncode != 0:
                    print("Error occurred while running the compiled executable.")
                    break

                it += 1
                with open('param.dat', 'r') as file:
                    lines = file.readlines()
                lines[-2] = f'myseed = 0{int(lines[-2].split("=")[-1])+1}\n'    # update random seed -> new random numbers   
                lines[-1] = f'iteration = {it}\n'                               # new interation value -> new folder 

                with open('param.dat', 'w') as file:
                    file.writelines(lines)

                os.remove(f'res/{keyword}/{box_x}_{delta_x}/{it-1}/V_{box_x}_pressure.dat')     # remove empty file

        run += 1    # -> change volume

    # change kind: -> from random to cubic
    with open('param.dat', 'r') as file:
        lines = file.readlines()
    lines[10] = f'kind = {ty}\n'
    with open('param.dat', 'w') as file:
        file.writelines(lines)

with open('param.dat', 'r') as file:
    lines = file.readlines()
for i in range(2, 5):
    lines[i] = f'box_{my_dict[i]} = {Volumes[0]}\n' 

lines[7] = f'delta_x = {all_delta_x[0]}\n'
lines[-2] = f'myseed = 02012024\n'  # reset seed    
lines[-1] = f'iteration = 1\n'      # reset counter
            
with open('param.dat', 'w') as file:
    file.writelines(lines)

