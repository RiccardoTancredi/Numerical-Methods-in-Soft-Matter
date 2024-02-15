import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy.optimize import curve_fit

plt.rcParams.update({'font.size': 13})  
folder = 'res/'
Temperatures = ['_0.5_T_c', '_T_c', '_2_T_c']
# 'T_c'
# types = ['energy', 'final_lattice_2D_Ising', 'magnetization']
# '.txt'

# make a color map of fixed colors
# cmap = colors.ListedColormap(['darkblue', 'darkred'])
# cmap = colors.LinearSegmentedColormap.from_list('my_colormap', ['darkblue', 'darkred'], 256)
cmap = colors.LinearSegmentedColormap.from_list('my_colormap', ['blue', 'red'], 256)

# Lattice plots
L = 50
lattice_plots = False

temp = [r'$0.5\cdot$', '', r'$2\cdot$']

if lattice_plots:
    for k, j in enumerate(Temperatures):
        for T in range(100, 1100, 100):
            file_name = folder+str(T)+j+'_lattice.txt'
            data = np.loadtxt(file_name)
            plt.figure(T)
            im = plt.pcolormesh(data, edgecolors='k', cmap=cmap, vmin=-1, vmax=1)
            plt.colorbar(im, cmap=cmap)
            plt.axis('off')
            name = temp[k] +r'$T_c$' + f' - Iteration = {T}'
            plt.title(name)
            ax = plt.gca()
            ax.set_aspect('equal')
            plt.savefig(f'img/lattice{j}_iteration_{T}_lattice.png', dpi=300)
            # plt.show(block=False)
            plt.clf()
    # plt.show()

cluster_size_hist = False

Temperatures = ['0.5_T_c', 'T_c', '2_T_c']
plt.clf()

if cluster_size_hist:
    for k, T in enumerate(Temperatures):
        fig, ax = plt.subplots(1)
        file_name = folder+str(T)+'_cluster_size.txt'
        data = np.loadtxt(file_name)
        n, bins, _ = ax.hist(data, label="data", density=True, bins=22, facecolor='mediumturquoise', edgecolor='black')
        ax.scatter((bins[1:]+bins[:-1])/2, n, marker='.', color='darkgreen')
        name = temp[k] +r'$T_c$'
        ax.set_title(f'{name}')
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$\mathcal{S}(\beta)$')
        ax.legend()
        ax.grid(ls='--', alpha=0.5)
        
        if k == 0:
            y_tot = [[0.00460, 0], [0.00460, 0], [0.00815, n[-1]], [0.00815, n[-1]]]
            x_tot = [[855, bins[-2]], [1920, bins[-1]], [1920, bins[-1]], [855, bins[-2]]]

            # Draw connecting lines from the corners of the last bin in the main plot to the corners of the inner plot figure
            for i in range(4):
                ax.plot(x_tot[i], y_tot[i], ls='--', color='0.8', lw=0.5)
            
            a = plt.axes([.4, .5, .3, .3])
            data_zoom = data[data >=(L**2-50)]
            n_zoom, bins_zoom, _ = plt.hist(data_zoom, label="last bin", density=True, bins=12, facecolor='mediumturquoise', edgecolor='black')
            plt.scatter((bins_zoom[1:]+bins_zoom[:-1])/2, n_zoom, marker='.', color='darkgreen')
            name = temp[k] +r'$T_c$' + f' - zoom in'
            plt.title(f'{name}')
            plt.xlabel(r'$x$')
            plt.ylabel(r'$\mathcal{S}(\beta)$')
            plt.legend()
            plt.grid(ls='--', alpha=0.5)
            a.set_xlim([min(data_zoom), max(data_zoom)])
            a.set_ylim([0, max(n_zoom) * 1.2])

        plt.savefig(f'img/{str(T)}_cluster_size.png', dpi=300)
        # plt.tight_layout()
        plt.show()

########################
# Autocorrelation time #
########################

# # t_eq:
# equilibrium_times = []
# # plt.clf()

# colors = ["#4B88A2", "#BB0A21"]
# labels=[r'$\mathcal{H}$', r'$m$']

# eq_times = True
# folder = 'res_autocorrelation/'
# for k, T in enumerate([20, 30, 50, 100, 150, 200]):
#     file_name_m = folder+f'{T}_T_c_magnetization.txt'
#     file_name_e = folder+f'{T}_T_c_energy.txt'
#     magn = (np.loadtxt(file_name_m)+1)/2
#     energy = np.loadtxt(file_name_e)
#     # Find the equilibium time
#     ind_e = np.where(abs(energy) >= abs((energy[::-1][:20]).mean()))[0][0]
#     ind_m = np.where(abs(magn) >= abs((magn[::-1][:20]).mean()))[0][0]
#     equilibrium_times.append(max([ind_e, ind_m]))
#     # print(equilibrium_times)

#     if eq_times:
#         plt.plot(energy[0:equilibrium_times[k]+1000], c=colors[0], label=labels[0])
#         plt.plot(magn[0:equilibrium_times[k]+1000], c=colors[1], label=labels[1])
#         # if T != '2_':
#         plt.vlines(equilibrium_times[k], energy[equilibrium_times[k]], 
#                 magn[equilibrium_times[k]], color='k', ls='--', label=r'$\tau_{eq}$')
            
#         plt.xlabel('MC steps per lattice site')
#         plt.ylabel('Observable')
#         # name = '.'.join(T.split("_"))[:-1] +r'$\cdot T_c$'
#         name = T
#         plt.title(name)
#         plt.legend()
#         plt.grid(ls='--', alpha=0.5)
#         # plt.savefig(f'img/energy_magn_{T}T_c.png', dpi=300)
#         plt.show()
#         # plt.clf()

# t_eq = equilibrium_times
# t_max = 10**3
# types = ['_magnetization', '_energy']
# correlation_function_e = np.zeros(t_max)
# correlation_function_m = np.zeros(t_max)
# corr_time_e, corr_time_m = [], []
# show_plots = False
# statistical_error_m = []

# for k, T in enumerate(Temperatures):
#     for j in range(2):
#         file_name = folder+T+types[j]+'.txt'
#         data = np.loadtxt(file_name) if j == 1 else np.abs(np.loadtxt(file_name))
#         average_observable = np.sum(data)/(len(data))
#         if j == 1:
#             # energy.append(average_observable)
#             for t in range(t_max):
#                 correlation_function_e[t] = np.sum(data[:t_max-t]*data[t:t_max])/(t_max-t)-((1./(t_max-t))**2)*np.sum(data[0:t_max-t])*np.sum(data[t:t_max])
#             correlation_function_e /= correlation_function_e[0]
#             # Now that we have our correlation function: -> estimate t_O = integrated autocorrelation time
#             # 2st method
#             tau = np.sum(correlation_function_e/correlation_function_e[0])
#             if t_max < tau:
#                 print('=================================================================')
#                 print('Correlations are present: re-evaluate the variance of the method!!')
#                 print('=================================================================')
#             corr_time_e.append(tau)
#             print(f'The integrated correlation time is tau = {tau}')
#             if show_plots:
#                 plt.plot(correlation_function_e, label='Correlation function')
#                 plt.show()
#             np.savetxt(f'auto_corr_e_{T}.txt', correlation_function_e)
#             # statistical_error_e.append(np.sqrt(np.sum((data[t_eq[k]:]-average_observable)**2)/(len(data)-t_eq[k]-1)))
#         else:
#             # magn.append(average_observable)
#             # X.append(np.var(data[t_eq[k]:])/(k_B*real_T)*L**2)
#             for t in range(t_max):
#                 correlation_function_m[t] = np.sum(data[0:t_max-t]*data[t:t_max])/(t_max-t)-(1./(t_max-t))**2 * np.sum(data[0:t_max-t])*np.sum(data[t:t_max])
#             correlation_function_m /= correlation_function_m[0]
#             # Now that we have our correlation function: -> estimate t_O = integrated autocorrelation time
#             # 2st method
#             tau = np.sum(correlation_function_m/correlation_function_m[0])
#             if t_max < tau:
#                 print('=================================================================')
#                 print('Correlations are present: re-evaluate the variance of the method!!')
#                 print('=================================================================')
#             print(f'The integrated correlation time is tau = {tau}')
#             corr_time_m.append(tau)
#             if show_plots:
#                 plt.plot(correlation_function_m, label='Correlation function')
#                 plt.show()
#             np.savetxt(f'auto_corr_m_{T}.txt', correlation_function_m)
#             statistical_error_m.append(np.sqrt(np.sum((data[t_eq[k]:]-average_observable)**2)/(len(data)-t_eq[k]-1)))

# n_uncor = [t_max/max(corr_time_e[i], corr_time_m[i]) for i in range(len(corr_time_m))]

# print(f'The correlation time for H = {corr_time_e}')
# print(f'The correlation time for m = {corr_time_m}')
# print(f'The number of uncorrelated measures = {n_uncor}')

# corr_time = [max(corr_time_e[i], corr_time_m[i]) for i in range(len(Temperatures))]


################
# z_w estimate #
################ 

estimate_z_W = True
corr_time_e, corr_time_m = [], []
show_plots = False
t_max = 10**3
types = ['_magnetization', '_energy']
correlation_function_e = np.zeros(t_max)
correlation_function_m = np.zeros(t_max)

print('=============')
print('z_W Estimate!')
print('=============\n')
                
if estimate_z_W:
    folder = 'res_autocorrelation/prova4/'
    taus = []
    LL = np.array([30, 50, 100, 150, 200])
    
    for k, L in enumerate(LL):
        for j in range(2):
            data = np.loadtxt(folder+f'{L}_cluster_size.txt')
            av = np.mean(data)
            mask = data < av
            av_cluster_size = np.mean(data[mask])
            file_name = folder+f'{L}'+types[j]+'.txt'
            data = np.loadtxt(file_name) 
            if j == 1:
                for t in range(t_max):
                    correlation_function_e[t] = np.sum(data[:t_max-t]*data[t:t_max])/(t_max-t)-((1./(t_max-t))**2)*np.sum(data[0:t_max-t])*np.sum(data[t:t_max])
                correlation_function_e /= correlation_function_e[0]
                # Now that we have our correlation function: -> estimate t_O = integrated autocorrelation time
                # 2st method
                tau = np.sum(correlation_function_e)
                corr_time_e.append(tau)
                print(f'The integrated correlation time for the energy for L = {L} is tau = {tau}')
                if show_plots:
                    plt.plot(correlation_function_e, label='Correlation function')
                    plt.legend()
                    plt.show()
            else:
                for t in range(t_max):
                    correlation_function_m[t] = np.sum(data[0:t_max-t]*data[t:t_max])/(t_max-t)-(1./(t_max-t))**2 * np.sum(data[0:t_max-t])*np.sum(data[t:t_max])
                correlation_function_m /= correlation_function_m[0]
                # Now that we have our correlation function: -> estimate t_O = integrated autocorrelation time
                # 2st method
                tau = np.sum(correlation_function_m)
                corr_time_m.append(tau)
                print(f'The integrated correlation time for the magnetization for L = {L} is tau = {tau}')
                if show_plots:
                    plt.plot(correlation_function_m, label='Correlation function')
                    plt.legend()
                    plt.show()
        taus.append(max(corr_time_m[k], corr_time_e[k])*av_cluster_size/L**2)
    
    print(f'taus = {taus}')
    popt, pcov = curve_fit(lambda x, a, b: a+x*b, np.log(LL), np.log(taus), p0=[1, 0.25])
    # I expect b == z_W
    print(f'popt = \n {popt}\n pcov = \n {pcov}\n')
    theor = (lambda x, a, b: np.exp(a)*x**b)(np.linspace(min(LL), max(LL)), *popt)
    

    # Metropolis part
    taus_Metropolis = []
    show_plots = False
    t_max = 10**4
    correlation_function_m = np.zeros(t_max)
    LL = np.array([30, 50, 100, 150, 200])
    
    folder_old = '../../L08/Exercise 3/finite_size_scaling/res/'
    for k, L in enumerate(LL):
        file_name = folder_old+f'{L}'+types[0]+'.txt'
        data = np.loadtxt(file_name)
        data = (data+1)/2
        t_max = 10**3 if L == 30 or L == 50 else 10**4
        correlation_function_m = np.zeros(t_max)
        for t in range(t_max):
            correlation_function_m[t] = np.sum(data[0:t_max-t]*data[t:t_max])/(t_max-t)-(1./(t_max-t))**2 * np.sum(data[0:t_max-t])*np.sum(data[t:t_max])
            
        correlation_function_m /= correlation_function_m[0]
        # Now that we have our correlation function: -> estimate t_O = integrated autocorrelation time
        if L != 150:
            # 2st method
            tau = np.sum(correlation_function_m)
        else:
            # 1st method -- better estimate in this case, but still an approximation
            popt_150, pcov_150 = curve_fit(lambda t, a, tau: a*np.exp(-t/tau), list(range(500)), correlation_function_m[:500])
            tau = popt_150[1]
        print(f'The integrated correlation time for the magnetization for L = {L} is tau = {tau}')
        if show_plots:
            plt.plot(correlation_function_m, label='Correlation function')
            plt.legend()
            plt.show()

        taus_Metropolis.append(tau)

    popt_Metr, pcov_Metr = curve_fit(lambda x, a, b: a+x*b, np.log(LL), np.log(taus_Metropolis), p0=[1, 2])
    # I expect b == z
    print(f'popt = \n {popt_Metr}\n pcov = \n {pcov_Metr}\n')
    theor_Metropolis = (lambda x, a, b: np.exp(a)*x**b)(np.linspace(min(LL), max(LL)), *popt_Metr)
    
    plt.scatter(LL, taus, label='Wolff', c='cornflowerblue')
    plt.plot(np.linspace(min(LL), max(LL)), theor, c='royalblue')
    plt.scatter(LL, taus_Metropolis, label='Metropolis', c='seagreen')
    plt.plot(np.linspace(min(LL), max(LL)), theor_Metropolis, c='forestgreen')
    
    plt.plot([], [], ' ', label=r'$z_W ='+f'{round(popt[1], 2)}'+r'$'+f'\n'+r'$z_M ='+f'{round(popt_Metr[1], 2)}'+r'$')
    
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(ls='--', alpha=0.5)
    plt.xlabel(r'$L$')
    plt.ylabel(r'$\tau$')
    plt.savefig(f'z_W_estimate.png', dpi=300)
    plt.show()