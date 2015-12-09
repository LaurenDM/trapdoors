import numpy as np
from ring_trapdoors import *
from gaussian import *
from joblib import Parallel, delayed

def preimage_sample_G_1D(k, u, s, n):
    x = []
    for i in range(k):
        if u % 2 == 0:
            x_i = gauss_samp_1D_even(s, u, n)
        else:
            x_i = gauss_samp_1D_odd(s, u, n)
        x.append(x_i)
        u = (u-x_i) / 2
    return np.array(x)


#   k: log q
#   u: target syndrome (ring element)
#   s: sigma
def process_input(k, u, s, i):
    return preimage_sample_G_1D(k, u[i], s, len(u))
def ring_preimage_sample_G(k, u, s):
    x = np.zeros((k, len(u)))
    parallel_output = Parallel(n_jobs=16)(delayed(process_input)(k, u, s, i) for i in range(len(u)))
    for i in range(len(u)):
        x[:,i] = parallel_output[i]
    return x

def ring_sample(sigma, mean, n):
    x = []
    for i in range(n):
        x.append(gauss_samp_1D(sigma, mean, n))
    return np.matrix(x)

def beta(x):
    return Rot(x)[0]

def get_rootSigma(A,R,E,s,q,r):
    (kplus2, n) = A.shape
    (k, n) = E.shape

    R1 = np.hstack([Rot(E[i]) for i in range(k)])
    R2 = np.hstack([Rot(R[i]) for i in range(k)])

    betaE = np.vstack([beta(E[i]) for i in range(k)])
    betaR = np.vstack([beta(R[i]) for i in range(k)])
    COV = r**2*np.vstack(
            [
            np.hstack(
                [Rot(np.matrix(A_mult(q, E, betaE))), Rot(np.matrix(A_mult(q, E, betaR))), R1]
                ),
            np.hstack(
                [Rot(np.matrix(A_mult(q, R, betaE))), Rot(np.matrix(A_mult(q, R, betaR))), R2],
                ),
            np.hstack(
                [R1.T, R2.T, np.eye(R1.shape[1])],
                ),
            ])

    SigmaP = s**2*np.eye(COV.shape[0]) - COV
    rootSigma = np.linalg.cholesky(SigmaP - (r/2)**2*np.eye(SigmaP.shape[0]))
    return rootSigma

def normal(mean, sigma):
    return np.random.normal(mean, sigma)

def preimage_sample_A(A,R,E,rootSigma, u, q, r):
    (kplus2, n) = A.shape
    k=kplus2-2

    floats = Parallel(n_jobs=16)(delayed(normal)(0,1) for i in range(kplus2 * n))
    #   perturbation: m = w+mbar ring elements
    p = rootSigma * np.matrix(floats).T
    parallel_output = Parallel(n_jobs=16)(delayed(gauss_samp_1D)(r/2, np.asscalar(x), n) for x in np.nditer(p))
    p = np.matrix(np.array(parallel_output))

    p = np.array(p.reshape(kplus2, n))

    v = A_mult(q, A, p)
    preimage = ring_preimage_sample_G(k, u-v, r)
    return p + combine_sample(q, R, E, preimage)
