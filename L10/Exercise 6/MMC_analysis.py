import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
plt.rcParams.update({'font.size': 13})

folder = 'MMC/'

T_c = 2./(np.log(1+np.sqrt(2))); 
T_vals = [T_c/2., 0.8*T_c, T_c, 1.5*T_c, 2.*T_c]
T_list = ['0.5_T_c_', '0.8_T_c_', 'T_c_', '1.5_T_c_', '2_T_c_']
T_latex = [r'$T = 0.5 \cdot T_c$', r'$T = 0.8 \cdot T_c$', 
           r'$T = T_c$', r'$T = 1.5 \cdot T_c$', r'$T = 2\cdot T_c$']
keys = ['lattice', 'magnetization', 'energy', 'swapping_rates']
ex = '.txt'
L = 50
cmap = colors.LinearSegmentedColormap.from_list('my_colormap', ['blue', 'red'], 256)
colors_0 = ['mediumseagreen', 'mediumturquoise', 'mediumblue', 
            'mediumslateblue', 'mediumpurple', 'mediumorchid']
colors = ['#067bc2', '#84bcda', '#ecc30b', '#f37748', '#d56062']   

all_energy = []

lattice_and_magn_ener_plots = True
if lattice_and_magn_ener_plots:
    # lattice plot
    for k, T in enumerate(T_vals):
        for c, it in enumerate(['', str(10**6)+'_']):
            file_name = folder+it+T_list[k]+keys[0]+ex
            data = np.loadtxt(file_name)
            plt.figure(k)
            im = plt.pcolormesh(data, edgecolors='k', cmap=cmap, vmin=-1, vmax=1)
            plt.colorbar(im, cmap=cmap)
            plt.axis('off')
            iter = r'$10^6$' if it == str(10**6)+'_' else r'$0$'
            name = f'L = {L} - '+ T_latex[k] + f' - Iteration = {iter}'
            plt.title(name)
            ax = plt.gca()
            ax.set_aspect('equal')
            plt.savefig(f'img/{folder}lattice{T_list[k]}_iteration_{it}_lattice.png', dpi=300)
            plt.show(block=True)
        plt.show()

        # magn and energy plot:
        file_name = folder+T_list[k]+keys[2]+ex
        data = np.loadtxt(file_name)
        all_energy.append(data)
        plt.plot(data[-1000:], label=r'$\mathcal{H}$', color='tab:blue')
        plt.hlines(np.mean(data), 0, 1000, ls='--', color='black')
        
        file_name = folder+T_list[k]+keys[1]+ex
        data = np.loadtxt(file_name)
        plt.plot(data[-1000:], label=r'$\mathcal{m}$', color='tab:red')
        plt.hlines(np.mean(data), 0, 1000, ls='--', color='black')
        plt.legend()
        plt.ylabel(f'Observable')
        plt.xlabel(f'MC steps per lattice size')
        plt.title('Magnetization & Energy time series - ' + T_latex[k])
        plt.tight_layout()
        plt.savefig(f'img/{folder}{T_list[k]}energy_magnetization.png', dpi=300)
        plt.show()

        
    # Energy distribution
    n_bins = 50
    for k, energy in enumerate(all_energy):
        plt.hist(energy, label=f'{T_latex[k]}', density=True, 
                bins=n_bins, facecolor=colors[k], 
                edgecolor='black', alpha=0.7, lw=1.2)

    plt.title(f'Energy distribution of adjacent chains')
    plt.legend()
    plt.grid(ls='--')
    plt.xlabel(r'$\mathcal{H}$')
    plt.ylabel(r'$\mathcal{P}\:(\mathcal{H})$')
    # plt.xlim(min(min(dataset) for dataset in all_energy), max(max(dataset) for dataset in all_energy))
    plt.tight_layout()
    plt.savefig(f'img/{folder}E_distr_chains.png', dpi=300)
    plt.show()


# Autocorrelation times
t_max = 10**3
correlation_function_e = np.zeros(t_max)
correlation_function_m = np.zeros(t_max)
corr_time_e, corr_time_m = [], []
show_plots = False

for k, T in enumerate(T_list):
    for j in range(1, 3):
        file_name = folder+T+keys[j]+ex
        data = np.loadtxt(file_name) # if j == 2 else np.abs(np.loadtxt(file_name))
        average_observable = np.sum(data)/(len(data))
        if j == 2:
            # energy.append(average_observable)
            for t in range(t_max):
                correlation_function_e[t] = np.sum(data[:t_max-t]*data[t:t_max])/(t_max-t)-((1./(t_max-t))**2)*np.sum(data[0:t_max-t])*np.sum(data[t:t_max])
            correlation_function_e /= correlation_function_e[0]
            # Now that we have our correlation function: -> estimate t_O = integrated autocorrelation time
            # 2st method
            tau = np.sum(correlation_function_e/correlation_function_e[0])
            corr_time_e.append(tau)
            print(f'The integrated correlation time is tau = {tau}')
            if show_plots:
                plt.plot(correlation_function_e, label='Correlation function')
                plt.legend()
                plt.show()
        else:
            # magn.append(average_observable)
            # X.append(np.var(data[t_eq[k]:])/(k_B*real_T)*L**2)
            for t in range(t_max):
                correlation_function_m[t] = np.sum(data[0:t_max-t]*data[t:t_max])/(t_max-t)-(1./(t_max-t))**2 * np.sum(data[0:t_max-t])*np.sum(data[t:t_max])
            correlation_function_m /= correlation_function_m[0]
            # Now that we have our correlation function: -> estimate t_O = integrated autocorrelation time
            # 2st method
            tau = np.sum(correlation_function_m/correlation_function_m[0])
            print(f'The integrated correlation time is tau = {tau}')
            corr_time_m.append(tau)
            if show_plots:
                plt.plot(correlation_function_m, label='Correlation function')
                plt.legend()
                plt.show()
            
n_uncor = [t_max/max(corr_time_e[i], corr_time_m[i]) for i in range(len(corr_time_m))]

print(f'The correlation time for H = {corr_time_e}')
print(f'The correlation time for m = {corr_time_m}')
print(f'The number of uncorrelated measures = {n_uncor}')

corr_time = [max(corr_time_e[i], corr_time_m[i]) for i in range(len(T_list))]

np.savetxt(f'Corr_time_swap.txt', np.array(corr_time))
np.savetxt(f'n_uncorr_swap.txt', np.array(n_uncor))

# Same for the 5 chains with swapping not allowed
correlation_function_e = np.zeros(t_max)
correlation_function_m = np.zeros(t_max)
corr_time_e, corr_time_m = [], []
show_plots = False
folder_no_swap = 'MMC_no_swap/'

for k, T in enumerate(T_list):
    for j in range(1, 3):
        file_name = folder_no_swap+T+keys[j]+ex
        data = np.loadtxt(file_name) # if j == 2 else np.abs(np.loadtxt(file_name))
        average_observable = np.sum(data)/(len(data))
        if j == 2:
            # energy.append(average_observable)
            for t in range(t_max):
                correlation_function_e[t] = np.sum(data[:t_max-t]*data[t:t_max])/(t_max-t)-((1./(t_max-t))**2)*np.sum(data[0:t_max-t])*np.sum(data[t:t_max])
            correlation_function_e /= correlation_function_e[0]
            # Now that we have our correlation function: -> estimate t_O = integrated autocorrelation time
            # 2st method
            tau = np.sum(correlation_function_e/correlation_function_e[0])
            corr_time_e.append(tau)
            print(f'The integrated correlation time is tau = {tau}')
            if show_plots:
                plt.plot(correlation_function_e, label='Correlation function')
                plt.legend()
                plt.show()
        else:
            # magn.append(average_observable)
            # X.append(np.var(data[t_eq[k]:])/(k_B*real_T)*L**2)
            for t in range(t_max):
                correlation_function_m[t] = np.sum(data[0:t_max-t]*data[t:t_max])/(t_max-t)-(1./(t_max-t))**2 * np.sum(data[0:t_max-t])*np.sum(data[t:t_max])
            correlation_function_m /= correlation_function_m[0]
            # Now that we have our correlation function: -> estimate t_O = integrated autocorrelation time
            # 2st method
            tau = np.sum(correlation_function_m/correlation_function_m[0])
            print(f'The integrated correlation time is tau = {tau}')
            corr_time_m.append(tau)
            if show_plots:
                plt.plot(correlation_function_m, label='Correlation function')
                plt.legend()
                plt.show()
            
n_uncor_no_swap = [t_max/max(corr_time_e[i], corr_time_m[i]) for i in range(len(corr_time_m))]

print(f'The correlation time for H (no swaps) = {corr_time_e}')
print(f'The correlation time for m (no swaps) = {corr_time_m}')
print(f'The number of uncorrelated measures (no swaps) = {n_uncor}')

corr_time_no_swap = [max(corr_time_e[i], corr_time_m[i]) for i in range(len(T_list))]

np.savetxt(f'Corr_time_NO_swap.txt', np.array(corr_time_no_swap))
np.savetxt(f'n_uncorr_NO_swap.txt', np.array(n_uncor_no_swap))

# Record swapping rates between chains
file_name = folder+keys[-1]+ex
data = np.loadtxt(file_name)
c_from = data[:, 0]
c_to = data[:, 1]
escape_rates = []
in_rates = []

for i in range(1, 6):
    # count the number of times from the i-th chain we go to the
    # j-th chain (not viceversa): I count only the 'escaping' rate, from
    # the i-th chain to any other one 
    i_mask = c_from == i
    escape_rate_from_i = sum(i_mask)/len(data)
    escape_rates.append(escape_rate_from_i)
    # count the number of times from a generic state j-th we go to the
    # i-th chain (not viceversa): I count only the 'ingoing' rate, from
    # a generic chain to the i-th one.  
    j_mask = c_to == i
    in_rates_to_j = sum(j_mask)/len(data)
    in_rates.append(in_rates_to_j)

#######################
# Swapping rates plot #
#######################

for i, el in enumerate(escape_rates):
    plt.scatter(i+1, el, label=r'$\lambda_'+f'{i+1}'+r'$', color=colors[i], marker='o')
plt.hlines(20/100, 1, 5, ls=':', color='k')
plt.title('Escape swapping rates per chain')
plt.xlabel('Chain number')
plt.ylabel(f'Escape rates ' + r'$\lambda$')
plt.xticks(list(range(1, 6)), list(range(1, 6)))
plt.grid(ls='--')
plt.legend()
plt.tight_layout()
plt.savefig(f'img/{folder}escape_rate_chain.png', dpi=300)
plt.show()

for i, el in enumerate(in_rates):
    plt.scatter(i+1, el, label=r'$p_{i'+f'{i+1}'+r'}$', color=colors[i], marker='o')
# plt.scatter(list(range(1, 6)), in_rates, label=r'$p_{ij}$', color='darkcyan', marker='o')
plt.hlines(20/100, 1, 5, ls=':', color='k')
plt.title(f'Ingoing swapping rates per chain '+r'$p_{ij}$')
plt.xlabel('Chain number')
plt.ylabel(f'Ingoing rates ' + r'$p_{ij}$')
plt.xticks(list(range(1, 6)), list(range(1, 6)))
plt.grid(ls='--')
plt.legend()
plt.tight_layout()
plt.savefig(f'img/{folder}in_rate_chain.png', dpi=300)
plt.show()

#########################
# Correlation time plot #
#########################

plt.scatter(T_vals, corr_time_no_swap, label=r'$\tau_{NO\: swap}$', 
            color=colors_0[0], marker='o')
plt.scatter(T_vals, corr_time, label=r'$\tau_{swap}$', 
            color=colors_0[-3], marker='o', s=80, alpha=0.5)
plt.title(f'Correlation time comparison')
plt.xlabel(r'$T_c$')
plt.ylabel(r'$\tau$')
T_list = ['0.5', '0.8', '1', '1.5', '2']
plt.xticks(T_vals, T_list)
plt.grid(ls='--')
plt.legend()
plt.tight_layout()
plt.savefig(f'img/{folder}tau_comparison.png', dpi=300)
plt.show()