from numpy import log2, ceil, floor
import numpy as np
import random
# import matplotlib.pyplot as plt
import gramschmidt
from Crypto.Util import *

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
    # plt.plot(x,y, 'go')
    # plt.show()

lookup_tables = {}

def get_lookup_table(sigma, mean, n):
    def rho(x):
        return np.exp(-(x - mean)**2 / (2.0 * sigma**2)) / (sigma*np.sqrt(2*np.pi))
    if (sigma, mean, n) not in lookup_tables:
        table = {}
        t = ceil(log2(n))
        for i in range(int(floor(mean-sigma*t)), int(ceil(mean+sigma*t + 1))):
            table[i] = rho(i)
        lookup_tables[(sigma, mean, n)] = table
    return lookup_tables[(sigma, mean, n)]

def gauss_samp_1D(sigma, mean, n):
    density_table = get_lookup_table(sigma, mean, n)
    sigma = float(sigma)
    mean = float(mean)
    t = ceil(log2(n))
    while True:

        #x = number.getRandomRange(int(floor(mean-sigma*t)), int(ceil(mean+sigma*t + 1)))
        x = np.random.randint(int(floor(mean-sigma*t)), int(ceil(mean+sigma*t + 1)))
        gaussian_density = density_table[x]
        coin_toss = random.uniform(0, 1)
        if coin_toss < gaussian_density:
            return x

def gauss_samp_1D_even(sigma, mean, n):
    while True:
        x = gauss_samp_1D(sigma, mean, n)
        if x % 2 == 0:
            return x

def gauss_samp_1D_odd(sigma, mean, n):
    while True:
        x = gauss_samp_1D(sigma, mean, n)
        if x % 2 == 1:
            return x

def gauss_samp(B, sigma, mean, n, q):
    Bg = gramschmidt.gs(B)
    current_center =  mean

    e = np.zeros(n)

    #print max([np.linalg.norm(np.array(Bg[:,i].T)[0]) for i in  range(n)])
    print B
    #print sigma
    # iterate i = n - 1 ... 0
    for i in range(0,n)[::-1]:

        # get the i^th gram schmidt vector
        bg_i = np.array(Bg[:,i].T)[0]

        cp_i = np.dot(current_center.T, bg_i) / np.dot(bg_i, bg_i)
        sp_i = sigma/np.linalg.norm(bg_i)

        #print "s_i, c_i: ", sp_i, cp_i
        z_i = gauss_samp_1D(sp_i, cp_i, n) % q

        # get the i^th basis vector
        b_i = np.array(B[:,i].T)[0]

        zb = z_i*b_i

        e = np.mod(np.add(e, zb), q)

        current_center = np.mod(np.subtract(current_center, zb), q)
    return e
