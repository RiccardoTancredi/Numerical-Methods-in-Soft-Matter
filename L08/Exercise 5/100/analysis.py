import numpy as np
import matplotlib.pyplot as plt

t_eq = np.loadtxt('t_eq.txt').astype(np.int64)
folder = 'res/'
Temperatures = ['0_1_', '0_9_', '2_']
# 'T_c'
types = ['energy', 'magnetization']
# '.txt'
real_t_max = 10**6
energy, magn = [], []
C_v, X = [], []
statistical_error_e = []
statistical_error_m = []
T_c = 2./(np.log(1+np.sqrt(2)))
k_B = 1.
t_max = 10_000
L = 100
correlation_function_e = np.zeros(t_max)
correlation_function_m = np.zeros(t_max)
corr_time_e, corr_time_m = [], []
show_plots = False

for k, T in enumerate(Temperatures):
    real_T = float('.'.join(T.split('_'))[:-1])*T_c
    for j in range(2):
        file_name = folder+T+'T_c_'+types[j]+'.txt'
        data = np.loadtxt(file_name) if j == 0 else np.abs(np.loadtxt(file_name))
        average_observable = np.sum(data[t_eq[k]:])/(len(data)-t_eq[k])
        if j == 0:
            energy.append(average_observable)
            C_v.append(np.var(data[t_eq[k]:])/(k_B*real_T**2)*L**2)
            for t in range(t_max):
                correlation_function_e[t] = np.sum(data[:t_max-t]*data[t:t_max])/(t_max-t)-((1./(t_max-t))**2)*np.sum(data[0:t_max-t])*np.sum(data[t:t_max])
            correlation_function_e /= correlation_function_e[0]
            # Now that we have our correlation function: -> estimate t_O = integrated autocorrelation time
            # 2st method
            tau = np.sum(correlation_function_e)
            if real_t_max < tau:
                print('=================================================================')
                print('Correlation are present: re-evaluate the variance of the method!!')
                print('=================================================================')
            corr_time_e.append(tau)
            print(f'The integrated correlation time is tau = {tau}')
            if show_plots:
                plt.plot(correlation_function_e, label='Correlation function')
                plt.show()
            np.savetxt(f'auto_corr_e_{k}.txt', correlation_function_e)
            statistical_error_e.append(np.sqrt((1+2*tau)*np.sum((data[t_eq[k]:]-average_observable)**2)/(len(data)-t_eq[k]-1)))
        else:
            magn.append(average_observable)
            X.append(np.var(data[t_eq[k]:])/(k_B*real_T)*L**2)
            for t in range(t_max):
                correlation_function_m[t] = np.sum(data[0:t_max-t]*data[t:t_max])/(t_max-t)-(1./(t_max-t))**2 * np.sum(data[0:t_max-t])*np.sum(data[t:t_max])
            correlation_function_m /= correlation_function_m[0]
            # Now that we have our correlation function: -> estimate t_O = integrated autocorrelation time
            # 2st method
            tau = np.sum(correlation_function_m)
            if real_t_max < tau:
                print('=================================================================')
                print('Correlation are present: re-evaluate the variance of the method!!')
                print('=================================================================')
            print(f'The integrated correlation time is tau = {tau}')
            corr_time_m.append(tau)
            if show_plots:
                plt.plot(correlation_function_m, label='Correlation function')
                plt.show()
            np.savetxt(f'auto_corr_m_{k}.txt', correlation_function_m)
            statistical_error_m.append(np.sqrt((1+2*tau)*np.sum((data[t_eq[k]:]-average_observable)**2)/(len(data)-t_eq[k]-1)))

n_uncor = [t_max/max(corr_time_e[i], corr_time_m[i]) for i in range(len(corr_time_m))]

print(f'Equilibrium time = {t_eq}\n')
print(f'Energy = {energy}\n')
print(f'Magn = {magn}\n')
print(f'C_v = {C_v}\n')
print(f'X = {X}\n')
print(f'Correlation times E = {corr_time_e}\n')
print(f'Correlation times M = {corr_time_m}\n')
print(f'# of measures uncorrelated = {n_uncor}\n')
print(f'Statistical error E = {statistical_error_e}\n')
print(f'Statistical error M = {statistical_error_m}\n')

with open('final_results.txt', 'w') as f:
    f.write('T:\t')
    for T in Temperatures:
        f.write(f'{T}*T_c\t')
    f.write('\n')
    f.write('t_eq:\t')
    for t in t_eq:
        f.write(f'{t}\t')
    f.write('\n')
    f.write('Energy:\t')
    for e in energy:
        f.write(f'{e}\t')
    f.write('\n')
    f.write('Stat error Energy:\t')
    for st_e in statistical_error_e:
        f.write(f'{st_e}\t')
    f.write('\n')
    f.write('Correlation time - Energy:\t')
    for corr_t in corr_time_e:
        f.write(f'{corr_t}\t')
    f.write('\n')
    f.write('Magnetization:\t')
    for m in magn:
        f.write(f'{m}\t')
    f.write('\n')
    f.write('Correlation time - Magn:\t')
    for corr_t in corr_time_m:
        f.write(f'{corr_t}\t')
    f.write('\n')
    f.write('Stat error Magn:\t')
    for st_m in statistical_error_m:
        f.write(f'{st_m}\t')
    f.write('\n')
    f.write('C_v:\t')
    for c_v in C_v:
        f.write(f'{c_v}\t')
    f.write('\n')
    f.write('X:\t')
    for x in X:
        f.write(f'{x}\t')
    f.write('\n')
    f.write('# uncorrelated measures n:\t')
    for un in n_uncor:
        f.write(f'{un}\t')
    f.write('\n')
