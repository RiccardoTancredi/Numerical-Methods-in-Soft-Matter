import numpy as np
import os
import subprocess


def generate_init_configuration(low=0, high=1, size=(1, 1), sigma=1):
    rows, cols = size
    data = np.random.uniform(low=low, high=high, size=(1, cols))
    # if rows/high**3 > 0.9:
    #     # to speed up the process for high density
    #     print('Creating the lattice...')
    #     coordinates = np.arange(rows)
    #     x, y, z = np.meshgrid(coordinates, coordinates, coordinates)
    #     # Stack the 3D arrays to get the final lattice
    #     data = np.vstack((x.flatten(), y.flatten(), z.flatten())).T
    # else:   
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
    N, N_strep, box_x, box_y, box_z, T, sigma, delta_x = data[:8]
    iteration = data[-1]
    T = 2 if T == 2 else 0.9

    # with open('param.dat', 'r') as file:
    #     lines = file.readlines()
    #     delta_x = float(lines[7].split('=')[-1])
    # delta_x = delta_x if delta_x != 1. else int(1) 

    return N, N_strep, box_x, box_y, box_z, T, sigma, delta_x, iteration, start_file_name


Temperatures = ['0.9', '2'] #  

for t in Temperatures:

    print(f'Temperature T = {t}')
    file_name = '../LJ_T2.dat' if t == '2' else '../LJ_T09.txt'
    prof_data = np.loadtxt(file_name)
    density_prof = prof_data[:, 0]
    pressure_prof = prof_data[:, 1]
    with open('param.dat', 'r') as file:
        lines = file.readlines()
    lines[0] = f'N = 200\n'
    lines[5] = f'temperature = {t}\n'
    lines[11] = f'myseed = 02012024\n' # reset seed
    with open('param.dat', 'w') as file:
        file.writelines(lines)

    N, N_strep, box_x, box_y, box_z, T, sigma, delta_x, iteration, start_file_name = read_param()
    file_path_x = f'{start_file_name}'

    Volumes = np.unique(np.round((N/density_prof)**(1/3)).astype(int))

    run = 0
    while run < len(Volumes):
        V = Volumes[run]
        print(f'\nBuffering V = {V}\n')
        
        os.makedirs(f'res/pressure/T_{t}/V_{V}/', exist_ok=True)
        
        it = 1
        while it < 11:
            print(f'\nBuffering iteration = {it}')
            
            # Modify lines in the 'params.dat' file:
            with open('param.dat', 'r') as file:
                lines = file.readlines()
            
            my_dict = {2:'x', 3:'y', 4:'z'}
            for i in range(2, 5):
                lines[i] = f'box_{my_dict[i]} = {V}\n' 
            lines[7] = f'delta_x = {V}\n' # update d_max

            lines[11] = f'myseed = 0{int(lines[11].split("=")[-1])+1}\n'
            with open('param.dat', 'w') as file:
                file.writelines(lines)

            N, N_strep, box_x, box_y, box_z, T, sigma, delta_x, iteration, start_file_name = read_param()
            
            # Generate initial (x, y, z) positions of the N particles
            # if t == '2':
            # init_x = np.random.uniform(low=0, high=box_x, size=(N, 3)) 
            # else:
                # init_x = generate_init_configuration(low=0, high=box_x, size=(N, 3), sigma=sigma)
            # np.savetxt(file_path_x, init_x)
            # print('\nFile created!')

            # Call C program
            executable_file = "mc.out"  

            # Run the compiled executable
            run_command = f"wsl ./{executable_file}"
            run_result = subprocess.run(run_command, shell=True)

            if run_result.returncode != 0:
                print("Error occurred while running the compiled executable.")
                break

            press = np.loadtxt(f'res/pressure/T_{T}/V_{V}_pressure.dat')
            en = np.loadtxt(f'res/pressure/T_{T}/V_{V}_energy.dat')
            if -1 < press[-1] < 11: 
                # This is just a check -> can be neglected
                np.savetxt(f'res/pressure/T_{T}/V_{V}/{it}_pressure.dat', press)
                np.savetxt(f'res/pressure/T_{T}/V_{V}/{it}_energy.dat', en)
                it += 1
        
        with open('param.dat', 'r') as file:
            lines = file.readlines()
        lines[11] = f'myseed = 02012024\n' # reset seed
        with open('param.dat', 'w') as file:
            file.writelines(lines)

        run += 1
            
with open('param.dat', 'r') as file:
    lines = file.readlines()
lines[0] = f'N = 100\n'
lines[5] = f'temperature = 1\n' # reset temperature
with open('param.dat', 'w') as file:
    file.writelines(lines)