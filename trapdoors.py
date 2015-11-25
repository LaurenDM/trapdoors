#!/usr/bin/python
import numpy as np
from Crypto.Util import *
from gaussian import gauss_samp_1D, test_1D_samp
import gramschmidt

# https://docs.google.com/document/d/1vqXmpGedS7U152H0F89yHwobM_jP47ipPs0vOiX6rMg/edit
def gen_trap(n,q,m):
    k = int(np.ceil(np.log2(q)))
    m0=m-n*k
    A0 = np.matrix(np.random.randint(0,q,(n,m0)))
    g = [int(2**i) for i in range(k)]
    G = np.kron(np.eye(n),g)
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

def trapdoor_G(n,q,k):
    Tg = np.diagflat(2*np.ones(k))+np.diagflat(-1*np.ones(k-1),-1)
    Tg[:,-1] = [int(b) for b in np.binary_repr(q,k)[::-1]]
    T = np.kron(np.eye(n),Tg)
    return T

def gen_basis(n, q, m, A, R):
    k = int(np.ceil(np.log2(q)))
    m0=m-n*k
    M1 = np.concatenate((np.concatenate((np.eye(m0),np.zeros((n*k,m0))),axis=0),R),axis=1)
    Tg = trapdoor_G(n,q,k)
    A0 = np.mod(-A[:,:m0],q)
    W = binary_decomp(A0,k)
    M2 = np.concatenate((np.concatenate((np.eye(m0),W),axis=0),np.concatenate((np.zeros((m0,n*k)),Tg),axis=0)),axis=1)
    return np.mod(M1*M2,q)

if __name__ == "__main__":
    n=32
    q=2053
    k=int(np.ceil(np.log2(q)))
    m=2*n*k
    A,R = gen_trap(n,q,m)
    B = gen_basis(n,q,m,A,R)
    print np.amax(np.amax(np.mod(A*B,q)))

    Bg = gramschmidt.gs(B)
    print "Max Gram-Schmidt: ", max([np.linalg.norm(np.array(Bg[:,i].T)[0]) for i in  range(m)])
