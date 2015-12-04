import numpy as np
from Crypto.Util import *
from gaussian import gauss_samp_1D, test_1D_samp
import gramschmidt

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

def A_mult(q,A,p):
    kp2,n = A.shape
    u=np.zeros(n)
    for i in range(kp2):
        u = u + poly_mult(q,np.array(A[i])[0],p[i])
    return np.mod(u,q)

def poly_mult(q,p1,p2):
    rotp1 = np.matrix([np.roll(p1,j) for j in range(n)])
    return np.mod(np.dot(p2,rotp1),q)

def gen_trap(n,q):
    k = int(np.ceil(np.log2(q)))
    a1 = np.array(np.random.randint(0,q,n))
    s=6.6 #smoothing parameter of integer lattice for eps=2^(-200) = 6.6
    gsamp = lambda i,j: gauss_samp_1D(s,0,n)
    vgsamp = np.vectorize(gsamp)
    R = np.matrix(np.fromfunction(vgsamp, (k,n), dtype=int))
    E = np.matrix(np.fromfunction(vgsamp, (k,n), dtype=int))
    g = [int(2**i) for i in range(k)]
    one = np.zeros(n)
    one[0]=1
    rota1=np.matrix([np.roll(a1,i) for i in range(n)])
    A1 = np.vstack([np.array(g[i]*one - (np.dot(R[i],rota1) +E[i])) for i in range(k)])
    one = np.zeros(n)
    one[0]=1
    A0 = np.vstack((np.matrix(one),np.matrix(a1)))
    A = np.vstack((A0,A1))
    return np.mod(A,q),R,E

def sample_g(n,k,u):
    X = np.zeros((k,n))
    for i in range(len(u)):
        X[:,i]=[int(b) for b in np.binary_repr(u[i],k)[::-1]]
    return X

def combine_sample(r,e,x):
    p0 = np.vstack(( A_mult(q,e,x), A_mult(q,r,x) ))
    return np.vstack((p0,x))

if __name__ == "__main__":
    # n=128
    # q=2048
    # k = int(np.ceil(np.log2(q)))
    # A,r,e=gen_trap(n,q)
    # u = np.array(np.random.randint(0,q,n))
    # x = sample_g(n,k,u)
    # z = combine_sample(r,e,x)
    # #p = combine_sample(r,e,np.hstack((np.ones((k,1)),np.zeros((k,n-1)))))
    # testresult = A_mult(q,A,z)
    # print np.mod(testresult-u,q)

    a=np.matrix([[1,2,3,4,5]])
    print Rot(a)
