from math import log, ceil, floor
import numpy as np
import random
import matplotlib.pyplot as plt
import gramschmidt

def test_1D_samp():
    d = {}
    for i in xrange(100000):
       x = gauss_samp_1D(8, 1, 512)
       if x not in d:
           d[x] = 1
       else:
           d[x] += 1
    l = []
    for k in d:
        l.append((k,d[k]))
    l.sort(key = lambda kv: kv[0])
    y = map(lambda kv: kv[1], l)
    x = map(lambda kv: kv[0], l)
    plt.plot(x,y, 'go')
    plt.show()

lookup_tables = {}

def get_lookup_table(sigma, mean, n):
    def rho(x):
        return np.exp(-(x - mean)**2 / (2.0 * sigma**2)) / (sigma*np.sqrt(2*np.pi))
    if (sigma, mean, n) not in lookup_tables:
        table = {}
        t = ceil(log(n, 2))
        for i in range(int(floor(mean-sigma*t)), int(ceil(mean+sigma*t + 1))):
            table[i] = rho(i)
        lookup_tables[(sigma, mean, n)] = table
    return lookup_tables[(sigma, mean, n)]

def gauss_samp_1D(sigma, mean, n):
    density_table = get_lookup_table(sigma, mean, n)
    sigma = float(sigma)
    mean = float(mean)
    t = ceil(log(n, 2))
    while True:
        x = np.random.randint(int(floor(mean-sigma*t)), int(ceil(mean+sigma*t + 1)))
        gaussian_density = density_table[x]
        coin_toss = random.uniform(0, 1)
        if coin_toss < gaussian_density:
            return x


def gauss_samp(B, sigma, mean, n):
    Bg = gs(B)

    current_center =  mean

    e = np.zeros()

    # iterate i = n - 1 ... 0
    for i in range(0,n)[::-1]:
        bg_i = Bg[:,i]
        cp_i = np.dot(current_center, bg_i) / np.dot(bg_i, bg_i)
        sp_i = sigma[i]/np.linalg.norm(bg_i)
        z_i = gauss_samp_1D(sp_i, cp_i, n)

        b_i = B[:,i]
        zb = np.mul(z_i, b_i)

        e = np.add(e, zb)

        current_center = current_center - zb

    return e

if __name__ == "__main__":
    pass
