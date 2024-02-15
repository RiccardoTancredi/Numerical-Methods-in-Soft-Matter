import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib.colors import Normalize


# file_name = "output_ex_2_2_1.txt"
# file_name = "output_ex_2_2_2.txt"
# file_name = 'additional_3.txt'
file_name = 'find_best_p.txt'
# file_name = "wrong.txt"
# file_name = "correct_disk_sampling.txt"
# file_name = "gaussian_sample.txt"
# file_name = "rejection_method.txt"

data = np.loadtxt(file_name)
# data = data[data <= 10]
r = data[:, 0]
theta = data[:, 1]
# x = data[:, 0]
# xi = data[:, 1]
# ind = np.where(x>0.5)[0]
# x = x[ind]
# xi = xi[ind]

x = r*np.cos(theta)
y = r*np.sin(theta)

def rho_1(x, n):
    c = n+1
    return c*x**n

def rho_2(x, c):
    n = 2
    return c*x**n

def add_1(x, mu):
    return mu*np.exp(-mu*x)

def add_2(x):
    return 2*x*np.exp(-x**2)

def add_3(x, a, b, n):
    return 1/(a+b*x)**n

def gaussian(x, y):
    return np.exp(-(x**2 + y**2)/2)/(2*np.pi)

def g_PDF(x, p):
    A = (2.*p)/(2*p**2 + 1)
    if 0 <= x <= p:
        return A
    return A/p * x * np.exp(p**2 - x**2)

def f_PDF(x):
    return np.sqrt(2/np.pi) * np.exp(-x**2)


n = 3
# c = 3/8
mu = 1.
a = 1
b = -(a**(1-n))/(1-n)
# p = 0.2# 1/np.sqrt(2) #0.5
# theta_t = np.linspace(0, 2*np.pi, 1000)
# x_t, y_t = circle = [np.cos(theta_t), np.sin(theta_t)]

# x_t = np.linspace(0, 10, 10_000)
# y_t = add_3(x_t, a, b, n)

# A = (2.*p)/(2*p**2 + 1)
# c = np.sqrt(2/np.pi)/A
# c = 1/A
# /A * np.sqrt(2/np.pi)*np.exp(-p**2)
# c = [np.sqrt(2/np.pi)/A if x[i] <= p else 1/A * np.sqrt(2/np.pi)*np.exp(-p**2) for i in range(len(x))]
# print(c, A, np.sqrt(2/np.pi)*np.exp(-p**2))
# y = np.array([g_PDF(x[i], p)*xi[i]*c for i in range(len(x))])


# X_t, Y_t = np.meshgrid(x, y)
# z = gaussian(x, y)
# heatmap, xedges, yedges = np.histogram2d(x, y, bins=1000, weights=z)
# norm = Normalize(vmin=np.min(z), vmax=np.max(z))

# # plt.imshow(heatmap.T, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], origin='lower', cmap='Greens')
# plt.imshow(heatmap.T, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], origin='lower', cmap='Greens', aspect='auto', norm=norm)
# colorbar = plt.colorbar(label=r'$\mathcal{N}(0,1)$', format='%.3f')
# colorbar.set_label(r'$\mathcal{N}(0,1)$', rotation=45, labelpad=15)
# plt.title('2D Heatmap of Gaussian sample')
# plt.xlabel('x')
# plt.ylabel('y')
# plt.xlim((-4,4))
# plt.ylim((-4,4))
# plt.grid(ls='--', alpha=0.5)
# plt.tight_layout()
# plt.savefig(f'gaussian_heatmap.png', dpi=300)
# plt.show()

# y = rho_2(x, c)

# inter = 5000
# z = gaussian(x, y)

# fig = plt.figure()
# ax = plt.axes(projection='3d')

# ax.contour3D(x_t, y_t, z_t, 30, cmap='gray')
# ax.scatter3D(x[:inter], y[:inter], z[:inter], c=z[:inter], cmap='Greens', label=f'Samples: N = {inter}')
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel(r'$\mathcal{N}(0, 1)$')
# ax.set_title('Gaussian sample')
# plt.savefig("gaussian_sampling.png", dpi=300)
# plt.show()

# inter = 1_000_000
# plt.scatter(x[:inter], y[:inter], label=f'Samples: N ='+r'$10^6$'+ f' of '+r'$g(x)$')
# plt.plot(x_t, f_PDF(x_t)*np.sqrt(np.pi/2), c = 'r', ls='--', label=r'$f(x)=\sqrt{\frac{2}{\pi}}e^{-x^2}$')
# c = [np.sqrt(2/np.pi)/A if x_t[i] <= p else p for i in range(len(x_t))]
# y_g = [g_PDF(x_t[i], p)*c for i in range(len(x_t))]
# plt.plot(x_t, y_g, c = 'g', ls='--', label=r'$c\cdot g(x)$')
# plt.scatter(x[:3000], y[:3000])
# plt.hist(np.abs(y)[::-1], alpha=0.9, density=True, bins=12)
# print(n, bins)
# plt.scatter((bins[1:]+bins[:-1])/2, n, marker='o', color='dimgrey')
# plt.hist(y[y>0], label="", density=True, bins=100)
# plt.plot(x_t, y_t, label='Theoretical distribution')
# n, bins, _ = plt.hist(data, label='Samples', density=True)
# n, bins, _ = plt.hist(data, label="Samples", density=True, bins=50, facecolor='mediumturquoise', edgecolor='black')
# plt.scatter((bins[1:]+bins[:-1])/2, n, marker='.', color='darkgreen')
# plt.title('Sampled data from ' +r'$\rho=\frac{1}{(a+bx)^n}$')
# plt.xlabel(r'$x$')
# plt.ylabel(r'$\rho(x)$')
# plt.legend()
# plt.grid(ls='--', alpha=0.5)
# plt.savefig("add_3.png", dpi=300)
# plt.show()
# plt.text(2.2, 0.7, r'$c \geq \sqrt{\frac{2}{\pi}}\cdot \frac{1}{A}$')
# plt.text(2.2, 0.5, r'$p = 0.2$')
# plt.ylabel(r'$\mathcal{N}(0, 1)$')

# find best p
effic = data/10**6
effic = np.insert(effic, 0, 0)
p = np.arange(0, 3.1, 0.1)
best_p = round(p[np.where(effic == max(effic))[0]][0], 1)
plt.plot(p, effic, c='steelblue', label=r'$p$', ls='dashed')
plt.hlines(max(effic), 0, max(p), color='royalblue', ls='-', label=r'$p_{best}=$'+f'{best_p}')
plt.legend(loc='right')
plt.grid(ls='--', alpha=0.5)
plt.title(r'$p$'+f' vs Rejection method efficiency')
plt.xlabel(r'$p$')
plt.ylabel('Efficiency')
# plt.savefig("wrong_r_distribution.png", dpi=300)
plt.savefig("best_p_estimate.png", dpi=300)
plt.show()

print(f'The best p is p = {best_p}')