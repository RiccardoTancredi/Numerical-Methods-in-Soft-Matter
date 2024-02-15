import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 13})  
folder = 'res/'
Omega = [int(10**2), int(10**3), int(10**4)]
omega = Omega[0]
file_name = folder + f'{omega}_brusselator.txt'
data = np.loadtxt(file_name, dtype=int)

if omega == 100:    
    X = data[:int(10**5), 0]
    Y = data[:int(10**5), 1]
else:    
    X = data[:, 0]
    Y = data[:, 1]

dt = 10**(-3)
tot_time = len(X)
t = np.arange(0, tot_time)

plt.plot(t*dt, X, label=r'$X$', color='royalblue')
plt.plot(t*dt, Y, label=r'$Y$', color='tab:red')
plt.xlabel(r'$t\:[s]$')
plt.ylabel(r'$\vec{\mathcal{C}}$', rotation=0)
plt.title(r'$\Omega = 10^2$' +' - ' r'$\mathcal{C}(0) = ' + f'{tuple(data[0])}' + r'$')
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
plt.legend()
plt.grid(ls='--')
plt.tight_layout()
plt.savefig(f'img/brusselator_dynamics_{int(omega)}.png', dpi=300)
plt.show()

# data_1 = np.loadtxt(file_name_1, dtype=int)
# prey_1 = data_1[:, 0]
# predators_1 = data_1[:, 1]

plt.plot(X, Y, label=r'$\mathcal{C}(0) = ' + f'{tuple(data[0])}' + r'$', 
         color='cornflowerblue')
# plt.plot(prey_1, predators_1, label=r'$\mathcal{C}(0) = ' + f'{tuple(data_1[0])}' + r'$', 
#          color='teal')
# plt.ylim((0, 1500))
plt.xlabel(r'$X$')
plt.ylabel(r'$Y$')
plt.title('Brusselator limit cycle - ' + r'$\Omega = 10^2$')
plt.legend()
plt.grid(ls='--')
plt.tight_layout()
plt.savefig(f'img/brusselator_{int(omega)}.png', dpi=300)
plt.show()