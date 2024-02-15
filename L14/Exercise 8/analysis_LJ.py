import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 13})  
show_inter_plot = False

temperatures = ['2', '0.9'] #  
for T in temperatures:
        
    folder = f'res/pressure/T_{T}'

    file_name = '../LJ_T09.txt' if T == '0.9' else '../LJ_T2.dat'
    prof_data = np.loadtxt(file_name)
    density_prof = prof_data[:, 0]
    pressure_prof = prof_data[:, 1]
    N = 200 
    print(f'N = {N}')
    Volumes = np.round((N/density_prof)**(1/3)).astype(int)

    density, pressure = [], []
    yerr = []

    for V in set(Volumes):
        data, data_std = [], []
        for it in range(1, 11):
            data_p = np.loadtxt(f'{folder}/V_{V}/{it}_pressure.dat')
            index = np.where(data_p < np.mean(data_p))[0][0]
            print(np.mean(data_p[index:]))
            data.append(np.mean(data_p[index:]))
            data_std.append(np.std(data_p[index:]))
            
            if show_inter_plot:
                # pressure plot
                plt.plot(data_p)
                plt.vlines(index, min(data_p), max(data_p), 
                        ls=':', color='k')
                plt.title(f'V = {V} - it = {it}')
                plt.show()

        data = np.array(data)
        data_std = np.array(data_std)
        # print(f'V = {V}: {data} -> {data.mean()}')
        pressure.append(data.mean())
        yerr.append(np.sqrt(np.sum(data_std**2))/len(data_std))
        density.append(N/V**3)

    plt.errorbar(density, pressure, yerr=yerr, color='royalblue',
                label='Simulated data\n'+r'$T^\star=$'+f'{T}', fmt='o',
                markersize=6, capsize=6)
    plt.scatter(density_prof, pressure_prof, color='tab:red', 
                label='State equation', marker='o')
    plt.xlabel(r'$\rho$')
    plt.ylabel(r'$P$')
    plt.grid(ls='--', alpha=0.5)
    plt.legend()
    # plt.savefig(f'../img/LJ_T{T}.png', dpi=300)
    plt.show()