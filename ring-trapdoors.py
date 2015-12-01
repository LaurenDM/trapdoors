import numpy as np
from Crypto.Util import *
from gaussian import gauss_samp_1D, test_1D_samp
import gramschmidt

def poly_mult(A,p):
    kp2,n = A.shape
    u=np.zeros(n)
    for i in range(kp2):
        rota = np.matrix([np.roll(np.array(A[i])[0],j) for j in range(n)])
        u = u + np.dot(p[i],rota)
    return u

def gen_trap(n,q):
    k = int(np.ceil(np.log2(q)))
    a1 = np.array(np.random.randint(0,q,n))
    s=6.6 #smoothing parameter of integer lattice for eps=2^(-200) = 6.6
    gsamp = lambda i,j: gauss_samp_1D(s,0,n)
    vgsamp = np.vectorize(gsamp)
    R = np.matrix(np.fromfunction(vgsamp, (k,n), dtype=int))
    E = np.matrix(np.fromfunction(vgsamp, (k,n), dtype=int))
    g = [int(2**i) for i in range(k)]
    rota1=np.matrix([np.roll(a1,i) for i in range(n)])
    A1 = np.vstack([np.array(g[i] - (np.dot(R[i],rota1) +E[i])) for i in range(k)])
    one = np.zeros(n)
    one[0]=1
    A0 = np.vstack((np.matrix(one),np.matrix(a1)))
    A = np.vstack((A0,A1))
    return A,R,E

if __name__ == "__main__":
    n=16
    q=1024
    k = int(np.ceil(np.log2(q)))
    A,r,e=gen_trap(n,q)
    print A
    print r
    print e
    ptest0 = np.vstack((e.sum(axis=0),r.sum(axis=0)))
    ptest1 = np.concatenate((np.ones((k,1)),np.zeros((k,n-1))),axis=1)
    ptest = np.vstack((ptest0,ptest1))
    print poly_mult(A,ptest)
