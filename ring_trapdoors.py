import numpy as np
from Crypto.Util import *
from gaussian import gauss_samp_1D, test_1D_samp
import gramschmidt
from ring_gaussamp import *

# General comment: each row in a matrix = coefficients of 1 polynomial
# Matrix = 'vector' of polynomials

def Rot(a): # a must be a np.matrix
    n = a.shape[1]
    rot=np.zeros((n,n))
    rot[:,0]=a
    for i in range(1,n):
        a=np.roll(a,1)
        a[0,0]=-a[0,0]
        rot[:,i]=a
    return rot

#   component-wise polynomial multiplication
def A_mult(q,A,p):
    kp2,n = A.shape
    u=np.zeros(n)
    for i in range(kp2):
        u = u + np.mod(np.dot(Rot(A[i]),p[i].T),q)
    return np.mod(u,q)

#   sum_i A[i] \cdot r
def R_mult(q, A, r):
    kp2,n = A.shape
    u=np.zeros(n)
    for i in range(kp2):
        u = u + np.mod(np.dot(Rot(A[i]),r.T),q)
    return np.mod(u,q)


def gen_trap_lwe(n,q):
    k = int(np.ceil(np.log2(q)))
    a1 = np.array(np.random.randint(0,q,n))
    s=4.7 #smoothing parameter of integer lattice for eps=2^(-200) = 6.6
    gsamp = lambda i,j: gauss_samp_1D(s,0,n)
    vgsamp = np.vectorize(gsamp)
    R = np.matrix(np.fromfunction(vgsamp, (k,n), dtype=int))
    E = np.matrix(np.fromfunction(vgsamp, (k,n), dtype=int))
    g = [int(2**i) for i in range(k)]
    one = np.zeros(n)
    one[0]=1
    rota1=Rot(np.matrix(a1))
    A1 = np.vstack([np.array(g[i]*one - (np.dot(rota1,R[i].T).T +E[i])) for i in range(k)])
    one = np.zeros(n)
    one[0]=1
    A0 = np.vstack((np.matrix(one),np.matrix(a1)))
    A = np.vstack((A0,A1))
    return np.mod(A,q),R,E

def gen_trap_sis(n, mbar, q):
    k = int(np.ceil(np.log2(q)))
    a = np.matrix(np.array(np.random.randint(0,q,n*mbar))).reshape(mbar, n)
    s=4.7 #smoothing parameter of integer lattice for eps=2^(-200) = 6.6
    gsamp = lambda i,j: gauss_samp_1D(s,0,n)
    vgsamp = np.vectorize(gsamp)
    R = np.array(np.fromfunction(vgsamp, (k*mbar,n), dtype=int))
    g = [int(2**i) for i in range(k)]
    one = np.zeros(n)
    one[0]=1
    A1 = np.vstack([np.array(g[i]*one - A_mult(q, a, R[i*mbar:(i+1)*mbar]).T) for i in range(k)])
    A = np.vstack((a,A1))
    return np.mod(A,q),R

def sample_g(n,k,u):
    X = np.zeros((k,n))
    for i in range(len(u)):
        X[:,i]=[int(b) for b in np.binary_repr(u[i],k)[::-1]]
    return X

def combine_sample(q, r,e,x):
    p0 = np.vstack(( A_mult(q,e,x), A_mult(q,r,x) ))
    return np.vstack((p0,x))

if __name__ == "__main__":
    n=512
    q=2**8
    k = int(np.ceil(np.log2(q)))
    mbar = int(np.ceil(np.log2(n)))
    A,r=gen_trap_sis(n, mbar, q)
    u = np.array(np.random.randint(0,q,n))
    rootSigma = get_rootSigma_sis(A,r,10000,q,2)
    preimage = preimage_sample_A_sis(A, r, rootSigma, u, q, 2)
    print np.mod(A_mult(q, A, preimage) - u, q)
