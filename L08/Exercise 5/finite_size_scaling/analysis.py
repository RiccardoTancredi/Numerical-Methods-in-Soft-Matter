import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

LL = np.array([10, 20, 35, 50, 70, 80, 100, 150, 200])
name =  '_magnetization.txt'

T_c = 2./(np.log(1+np.sqrt(2)))
k_B = 1.
t_max = 1_000
folder = 'res/'

show_plots = True
M_over_V, X_over_V = [], []
fig, ax = plt.subplots(1, 3, figsize=(12, 3))
for k, L in enumerate(LL):
    file_name_m = folder + str(L) + name
    magn = np.loadtxt(file_name_m)
    # Find M(T) -> beta exponent
    
    M_over_V.append(np.abs(magn).mean())
    X_over_V.append(np.abs(magn).var()/(k_B*T_c))

    if show_plots and k < 3:
        ax[k].plot(magn[:t_max], color='lightseagreen', label="$m$")
        ax[k].set_title(r'$L='+f'{L}'+r',\: T=T_c$') 
        ax[k].set_xlabel('MC steps per lattice site')
        if k == 0:
            ax[k].set_ylabel(r'Observable')
        ax[k].hlines(M_over_V[k], 0, t_max, color='gold', label=r'$\langle m \rangle$')
        ax[k].grid(ls='--', alpha=0.5)
        ax[k].legend()
        # plt.plot(list(range(t_max)), magn[:t_max])
        # plt.show()

    print(f'The equilibrium M for L = {L} is {M_over_V[k]}')
    print(f'The equilibrium X for L = {L} is {X_over_V[k]}')
plt.tight_layout()
plt.savefig('img/magnetization_L.png', dpi=300)
plt.show()

############
# beta fit #
############
y = np.array(M_over_V)
popt, pcov = curve_fit(lambda x, a, b: a*x+b, np.log(LL), np.log(y), p0=[-1/8, 1])

print(f'The extimated popt = \n {popt} \n and pcov = \n {pcov}\n')

theor = (lambda x, a, b: b*x**a)(LL, popt[0], np.exp(popt[1]))
print(f'y = {y}')
print(f'Theor = {theor}')
plt.scatter(LL, y, label='data', c='cornflowerblue')
plt.plot(LL, theor, label='fit', c='darkred')
plt.plot([], [], ' ', label=r'$-\beta/\nu ='+f'{round(popt[0], 3)}'+r'$')
plt.legend()
plt.grid(ls='--', alpha=0.5)
plt.yscale('log')
plt.xscale('log')
plt.title(r'$\beta\:$ estimate', fontsize=16)
plt.xlabel(r'$L$', fontsize=18)
plt.ylabel(r'$\langle m \rangle$', fontsize=18)
plt.tight_layout()
plt.savefig('img/beta_estimate.png', dpi=300)
plt.show()

print(f'The expected -beta/nu = {-1/8} and here we get {popt[0]}')

with open('exponents_fitting.txt', 'w') as f:
    f.write('############\n# beta fit #\n############\n')
    f.write('a,\t b\n')
    for el in popt:
        f.write(f'{el}\t')
    f.write('\nCov matrix:\n')
    for i, el in enumerate(pcov):
        if i % 2 == 1:
            f.write('\n')
        f.write(f'{el}\t')



#############
# gamma fit #
#############
y = np.array(X_over_V*LL**2)
popt, pcov = curve_fit(lambda x, a, b: a*x+b, np.log(LL), np.log(y), p0=[7/4, 1])

print(f'The extimated popt = \n {popt} \n and pcov = \n {pcov}\n')

theor = (lambda x, a, b: b*x**a)(LL, popt[0], np.exp(popt[1]))
print(f'y = {y}')
print(f'Theor = {theor}')
plt.scatter(LL, y, label='data', c='cornflowerblue')
plt.plot(LL, theor, label='fit', c='darkred')
plt.plot([], [], ' ', label=r'$\gamma/\nu ='+f'{round(popt[0], 2)}'+r'$')
plt.legend()
plt.grid(ls='--', alpha=0.5)
plt.yscale('log')
plt.xscale('log')
plt.title(r'$\gamma\:$ estimate', fontsize=16)
plt.xlabel(r'$L$', fontsize=18)
plt.ylabel(r'$N\cdot\chi_T$', fontsize=18)
plt.tight_layout()
plt.savefig('img/gamma_estimate.png', dpi=300)
plt.show()

print(f'The expected gamma/nu = {7/4} and here we get {popt[0]}')

with open('exponents_fitting.txt', 'a') as f:
    f.write('############\n# gamma fit #\n############\n')
    f.write('c,\t d\n')
    for el in popt:
        f.write(f'{el}\t')
    f.write('\nCov matrix:\n')
    for i, el in enumerate(pcov):
        if i % 2 == 1:
            f.write('\n')
        f.write(f'{el}\t')