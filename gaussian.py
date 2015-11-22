from math import log, ceil
import numpy as np
import random
import matplotlib.pyplot as plt
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

def gauss_samp_1D(sigma, mean, n):
    sigma = float(sigma)
    mean = float(mean)
    t = ceil(log(n, 2))
    def rho(x):
        return np.exp(-(x - mean)**2 / (2.0 * sigma**2)) / (sigma*np.sqrt(2*np.pi))
    while True:
        x = np.random.randint(mean-sigma*t, mean+sigma*t + 1)
        gaussian_density = rho(x)
        coin_toss = random.uniform(0, 1)
        if coin_toss < gaussian_density:
            return x
