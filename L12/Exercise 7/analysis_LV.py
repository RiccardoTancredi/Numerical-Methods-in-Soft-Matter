import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy.optimize import curve_fit

plt.rcParams.update({'font.size': 13})  
folder = 'res/'
file_name_1 = folder + 'population_stat.txt'
file_name = folder + 'population_no_stat.txt'
data = np.loadtxt(file_name, dtype=int)
prey = data[:, 0]
predators = data[:, 1]

dt = 10**(-4)
tot_time = len(prey)
t = np.arange(0, tot_time)

# plt.plot(t*dt, prey, label=r'$X_1:$ ' +f'prey', color='royalblue')
# plt.plot(t*dt, predators, label=r'$X_2:$ ' +f'predators', color='tab:red')
# plt.xlabel(r'$t\:[s]$')
# plt.ylabel(r'$\vec{\mathcal{C}}$', rotation=0)
# plt.title(r'$\mathcal{C}(0) = ' + f'{tuple(data[0])}' + r'$')
# plt.ylim((0, 1400))
# if file_name == 'population_no_stat.txt':
#     wh = min(np.where(np.isclose(predators, 0, rtol=10, atol=10) == True)[0])
#     start = np.array([9, 600])
#     end = np.array([(t*dt)[wh], 120-predators[wh]])
#     plt.arrow(*start, *(end-start), fc='teal', ec='teal', 
#               lw=2., head_length=80, head_width=.2)
# else:
#     plt.hlines(np.mean(prey), 0, 10, color='royalblue', ls=':')
#     plt.hlines(np.mean(predators), 0, 10, color='tab:red', ls=':')
#     plt.text(9, np.mean(prey)+25, r'$\mu = $'+ f'{round(np.mean(prey))}', 
#              color='royalblue')
#     plt.text(9.1, np.mean(predators)+25, r'$\mu = $' + f'{round(np.mean(predators))}', 
#              color='tab:red')
#     print(f'Prey mean = {np.mean(prey)} \n Predators mean = {np.mean(predators)}')
# plt.legend()
# plt.grid(ls='--')
# plt.tight_layout()
# # plt.savefig(f'img/population_far_from_stationarity.png', dpi=300)
# # plt.savefig(f'img/population_near_stationarity.png', dpi=300)
# plt.show()

data_1 = np.loadtxt(file_name_1, dtype=int)
prey_1 = data_1[:, 0]
predators_1 = data_1[:, 1]

plt.plot(prey, predators, label=r'$\mathcal{C}(0) = ' + f'{tuple(data[0])}' + r'$', 
         color='cornflowerblue')
plt.plot(prey_1, predators_1, label=r'$\mathcal{C}(0) = ' + f'{tuple(data_1[0])}' + r'$', 
         color='teal')
plt.ylim((0, 1500))
plt.xlabel('Prey')
plt.ylabel('Predators')
plt.title('Lotka - Volterra')
plt.legend()
plt.grid(ls='--')
plt.tight_layout()
# plt.savefig(f'img/prey_predators_dynamics.png', dpi=300)
plt.show()