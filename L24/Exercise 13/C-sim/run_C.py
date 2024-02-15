import numpy as np
import os
import subprocess

param_filename = 'parameters.dat'
save_folder = 'res'
forces = list(range(0, 10))
for i, f_active in enumerate(forces):
    print(f'\nF = {f_active}\n')
    
    # Update lines in the 'params.dat' file:
    with open(param_filename, 'r') as file:
        lines = file.readlines()
    lines[18] = f'f_active = {f_active}\n'
    with open(param_filename, 'w') as file:
        file.writelines(lines)

    # Call C program
    executable_file = "G-2"  

    # Run the compiled executable
    run_command = f"wsl ./{executable_file} {param_filename}"
    run_result = subprocess.run(run_command, shell=True)

    if run_result.returncode != 0:
        print("Error occurred while running the compiled executable.")
        break

    # After the file creation, move it to another folder and change its name
    filename = 'CONF/' + os.listdir('CONF/')[0]
    traj = np.loadtxt(filename)
    np.savetxt(f'{save_folder}/traj_f_active_iter_{i}.dat', traj)
    filename = 'DATA/' + os.listdir('DATA/')[0]
    filename = 'DATA/' + os.listdir('DATA/')[0] + '/' + os.listdir(filename)[0]
    x_rel = np.loadtxt(filename)
    np.savetxt(f'{save_folder}/x-rel_f_active_iter_{i}.dat', x_rel)

    print('\nFile saved in proper folder\n')
    
with open(param_filename, 'r') as file:
    lines = file.readlines()
# reset initial values
lines[18] = f'f_active = {0}\n'
with open(param_filename, 'w') as file:
    file.writelines(lines)
        