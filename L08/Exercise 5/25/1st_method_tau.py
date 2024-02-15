import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

'''
The integrated correlation time is tau = 13.750939417214848
The integrated correlation time is tau = 77.17423664981143
The integrated correlation time is tau = 1.3892447799450174
'''

i = 0
file_name = f'auto_corr_e_{i}.txt'
data = np.loadtxt(file_name)
plt.plot(data)
plt.show()
data = data[data>0]
popt, pcov = curve_fit(lambda t, a, tau: a*np.exp(-t/tau), list(range(len(data))), data)

print(f'params = {popt},\n cov matrix = \n{pcov}')