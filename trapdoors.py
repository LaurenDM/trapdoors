#!/usr/bin/python
import numpy as np
from Crypto.Util import *
from gaussian import gauss_samp_1D, test_1D_samp


if __name__ == "__main__":
    pass


def gen_trap(n,q,m):
    k = int(np.ceil(np.log2(q)))
    m0=m-n*k
    A0 = np.matrix(np.random.randint(0,q,(n,m0)))
    g = [int(2**i) for i in range(k)]
    G = np.zeros((n,n*k), dtype=int)
    for i in range(n):
        G[i,i*k:(i+1)*k]=g
    s=4.7 #smoothing parameter of integer lattice for eps=2^(-100) = 4.7
    gsamp = lambda i,j: gauss_samp_1D(s,0,n)
    vgsamp = np.vectorize(gsamp)
    R = np.matrix(np.fromfunction(vgsamp, (m0,n*k), dtype=int))
    A = np.concatenate((A0,G-A0*R),axis=1)
    R = np.concatenate((R,np.eye(n*k)),axis=0)
    return A,R

gen_trap(128,2053,3000)

