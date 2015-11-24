#!/usr/bin/python
import numpy as np
from Crypto.Util import *
from gaussian import gauss_samp_1D, test_1D_samp

# https://docs.google.com/document/d/1vqXmpGedS7U152H0F89yHwobM_jP47ipPs0vOiX6rMg/edit
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
    A = np.concatenate((A0,np.mod(G-A0*R,q)),axis=1)
    R = np.concatenate((R,np.eye(n*k)),axis=0)
    return A,R

def binary_decomp(A0,k):
    m0 = A0.shape[1]
    W=np.zeros((1,m0))
    for i in range(len(A0)):
        nextrow=np.zeros((k,m0))
        for j in range(m0):
            nextrow[:,j] = [int(b) for b in np.binary_repr(A0[i,j],k)[::-1]]
        W = np.concatenate((W,nextrow),axis=0)
    return W[1:]

def trapdoor_G(n,k):
    Tg = np.diagflat(2*np.ones(k))+np.diagflat(-1*np.ones(k-1),-1)
    T = np.zeros((n*k,n*k), dtype=int)
    for i in range(n):
        T[i*k:(i+1)*k,i*k:(i+1)*k]=Tg
    return T

def gen_basis(n, q, m, A, R):
    k = int(np.ceil(np.log2(q)))
    m0=m-n*k
    M1 = np.concatenate((R,np.concatenate((np.eye(m0),np.zeros((n*k,m0))),axis=0)),axis=1)
    Tg = trapdoor_G(n,k)
    A0 = np.mod(-A[:,:m0],q)
    W = binary_decomp(A0,k)
    M2 = np.concatenate((np.concatenate((Tg,np.zeros((m0,n*k))),axis=0),np.concatenate((W,np.eye(m0)),axis=0)),axis=1)
    return np.mod(M1*M2,q)

if __name__ == "__main__":
    A,R = gen_trap(32,2053,1024)
    B = gen_basis(32,2053,1024,A,R)
    print np.mod(A*B,2053)
