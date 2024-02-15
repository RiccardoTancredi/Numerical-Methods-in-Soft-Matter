#############################
# Multiple Histogram Method #
#############################

import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess

plt.rcParams.update({'font.size': 13})  
folder = 'res/MMC/50/'
T_c = 2./(np.log(1+np.sqrt(2)))
eps = 0.01
T_vals = np.array([T_c-2.*eps, T_c-eps, T_c, T_c+eps, T_c+2.*eps])
T_list = ["T_c_2eps_", "T_c_eps_", "T_c_", "T_c_p_eps_", "T_c_p_2eps_"]
T_latex = [r'$T_c - 2\varepsilon$', r'$T_c-\varepsilon$', 
           r'$T_c$', r'$T_c + \varepsilon$', 
           r'$T_c + 2\varepsilon$']

# folder = "../../L10/Exercise 4/MMC_no_swap/"
# T_c = 2./np.log(1+np.sqrt(2))
# T_vals = np.array([0.5*T_c, 0.8*T_c, T_c, 1.5*T_c, 2.*T_c])
# T_list = ["0.5_T_c_", "0.8_T_c_", "T_c_", "1.5_T_c_", "2_T_c_"]
# T_latex = [r'$0.5\cdot T_c$', r'$0.8\cdot T_c$', 
#            r'$T_c$', r'$1.5\cdot T_c $', 
#            r'$2\cdot T_c $']


keys = ['lattice', 'magnetization', 'energy', 'swapping_rates']
ex = '.txt'
colors = np.array(["#ddfff7", "#93e1d8", "#ffa69e", "#aa4465", "#460b2e"])
colors = colors[::-1]

betas = 1./T_vals
L = 50
N = L**2

file_n =f'res/mhm/{L}/direct_average.txt'
average_energy = np.loadtxt(file_n)

# betas_range = np.linspace(min(betas), max(betas), num=num)
betas_range = np.loadtxt(f'res/mhm/{L}/betas.txt')

    
file_n =f'res/mhm/{L}/Z_k.txt'

if not os.path.isfile(file_n):
    cpp_file = "exercise10_MHM.cpp"
    extra_cpp_file = "save_data.cpp"
    executable_file = "Ex10_MHM"  

    # Compile C++ file using g++ (wsl before since I work on Windows)
    compile_command = f"wsl g++ {cpp_file} {extra_cpp_file} -o {executable_file}"
    compile_result = subprocess.run(compile_command, shell=True)

    if compile_result.returncode == 0:
        # Run the compiled executable
        run_command = f"wsl ./{executable_file}"
        run_result = subprocess.run(run_command, shell=True)

        if run_result.returncode != 0:
            print("Error occurred while running the compiled executable.")
    else:
        print("Error occurred while compiling the C++ code.")

Z_k = np.loadtxt(file_n)
beta_j = 1./T_vals

betas = np.loadtxt(f'res/mhm/{L}/betas.txt')
Z_beta = np.loadtxt(f'res/mhm/{L}/Z_beta.txt')
U_beta = np.loadtxt(f'res/mhm/{L}/U_beta.txt')

skip = True
if not skip:
    plt.scatter(1./betas, U_beta/N, color='royalblue', marker='o',
                label=r'$U(T)^{MHM}$', s=40, alpha=0.8)
    plt.scatter(T_vals, average_energy/N, color='tab:red', s=100,
                label=r'$\langle U(T_k) \rangle$', marker='*')
    plt.vlines(T_c, min(U_beta/N), max(U_beta/N), ls=':', 
               color='k', label=r'$T_c$')
    plt.xlabel(r'$T$')
    plt.ylabel(r'$U(T)$')
    # plt.xticks(T_vals, T_latex)
    plt.legend(frameon=False)
    plt.grid(ls='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(f'img/MHM.png', dpi=300)
    plt.show()


# Specific heat
L = np.array([20, 30, 50])
N = L**2

expected_C = []
for i, l in enumerate(L):
    C = []
    for k, T in enumerate(T_list):
        data = np.loadtxt(f'res/MMC/{l}/{T}{keys[2]}{ex}')*N[i]
        C.append(np.var(data)/(N[i]*(beta_j[k]**2)))
    
    expected_C.append(np.array(C))

expected_C = np.array(expected_C)

minimum, maximum = 100, 100

want_zoom = True
zoom = 'zoom/' if want_zoom else ''

for i, l in enumerate(L):
    betas = np.loadtxt(f'res/mhm/{l}/{zoom}betas.txt')
    Z_beta = np.loadtxt(f'res/mhm/{l}/{zoom}Z_beta.txt')
    U_beta = np.loadtxt(f'res/mhm/{l}/{zoom}U_beta.txt')
    C_beta = np.loadtxt(f'res/mhm/{l}/{zoom}C_beta.txt')
    lab_1 = r'$C(T)^{MHM}$' if i == 0 else ''
    lab_2 = r'$C(T_i)$' if i == 0 else ''
    minimum = min(C_beta) if min(C_beta) < minimum else minimum
    maximum = max(C_beta) if max(C_beta) > maximum else maximum
    # col = colors[i] if not want_zoom else 'teal'
    plt.plot(1./betas, C_beta, color=colors[i], ls='-.',
                label=lab_1, alpha=1)
    if not want_zoom:
        plt.text(1./betas[np.where(C_beta==max(C_beta))[0][0]]+0.01, 
                 max(C_beta)+5, f'L={l}', color=colors[i])
    else:    
        plt.scatter(T_vals, expected_C[i], color='tab:red', s=100,
                    label=lab_2, marker='*')
if not want_zoom:
    plt.vlines(T_c, minimum, maximum+10, ls=':', 
               color='k', label=r'$T_c$')
if want_zoom:
    plt.text(T_vals[1], 55, r'L=50')
    plt.text(np.mean([T_vals[2], T_vals[3]]), 50, r'L=30')
    plt.text(T_vals[-2], 45, r'L=20')
plt.xlabel(r'$T$')
plt.ylabel(r'$C(T)$')
if want_zoom:
    plt.xticks(T_vals, T_latex)
loc = 'upper left' if want_zoom else 'upper right'
plt.legend(frameon=False, loc=loc)
plt.grid(ls='--', alpha=0.5)
plt.tight_layout()
appendix = '_zoom' if want_zoom else ''
plt.savefig(f'img/specific_heat_MHM{appendix}.png', dpi=300)
plt.show()

