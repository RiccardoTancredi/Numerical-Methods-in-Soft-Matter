import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# file_name = "rejection_method.txt"

# data = np.loadtxt(file_name)

# def sigmaf_over_f(T):
#     return 1- np.exp(-T)

def sigmaf(T):
    return np.sqrt(np.exp(-T)*(1- np.exp(-T)))

def sigma_a_star(T, a):
    return np.sqrt(np.exp(-T*(2-a))/(a*(2-a)) - np.exp(-2*T))

def a_star(T):
    # print(T, ((T+1)-np.sqrt(T**2+1))/T)
    return min(1, ((T+1)-np.sqrt(T**2+1))/T)

def sigmaF(T):
    best_a = a_star(T)
    return sigma_a_star(T, best_a)

def sigmaf_over_sigmag(T):
    sigma_up = sigmaf(T)
    # best_a = a_star(T)
    sigma_down = sigmaF(T)
    return sigma_up/sigma_down

T = [3, 5, 10, 20]
first = [sigmaf(t)/np.exp(-t) for t in T]
second = [sigmaF(t)/np.exp(-t) for t in T]
third = [sigmaf_over_sigmag(t) for t in T]

df = pd.DataFrame(columns=["T", "mf", "sf", "rapp"])
df.T = T
df.mf = first
df.sf = second
df.rapp = third

print(df)

plt.plot(T, first, marker='o', label=r'$\frac{\sigma(f)}{\langle f \rangle}$')
plt.plot(T, second, marker='o', label=r'$\frac{\sigma(a^\star, F(x))}{\langle f \rangle}$')
plt.plot(T, third, marker='o', label=r'$\frac{\sigma(f)}{\sigma(a^\star, F(x))}$')
# plt.hlines(1, min(T), max(T), colors='r', ls='--', label=r'$\frac{\sigma(f)}{\sigma(a^\star, F(x))}$') 
font = {'size':15}
plt.title('Importance sampling', fontdict=font)
plt.xlabel('T', fontdict=font)
plt.ylabel(r'$f(x)$', fontdict=font)
plt.xticks(np.linspace(3, 21, 4), T)
# plt.text(2.2, 0.7, r'$c \geq \sqrt{\frac{2}{\pi}}\cdot \frac{1}{A}$')
# plt.text(2.2, 0.5, r'$p = 0.2$')
# plt.ylabel(r'$\mathcal{N}(0, 1)$')
plt.legend(prop=font)
plt.grid(alpha=0.5, ls='--')
plt.yscale('log')
# plt.savefig("wrong_uniform.png", dpi=300)
plt.savefig("importance_sampling.png", dpi=300)
plt.show()