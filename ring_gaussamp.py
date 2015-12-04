import numpy as np
from ring_trapdoors import *
def preimage_sample_G_1D(k, u, s):
    x = []
    for i in range(k):
        if u % 2 == 0:
            x_i = gauss_samp_1D_even()
        else:
            x_i = gauss_samp_1D_odd()
        x.append(x_i)
        u = (u-x_i) / 2
    return np.matrix(x)


#   k: log q
#   u: target syndrome (ring element)
#   s: sigma
def ring_preimage_sample_G(k, u, s):
    x = []
    for i in range(len(u)):
        x.append(preimage_sample_G_1D(k, u[i], s))
    return x

def ring_sample(sigma, mean, n):
    x = []
    for i in range(n):
        x.append(gauss_samp_1D(sigma, mean, n))
    return np.matrix(x)

def beta(x):
    return Rot(x)[0]

def preimage_sample_A(A, R, E, u, s, q):
    print R.shape
    print E.shape
    Sigma = s**2
    SigmaG = 2
    SigmaP = Sigma - R * SigmaG * R.T

    R1 = np.hstack([Rot(E[i]) for i in range(E.shape[0])])
    R2 = np.hstack([Rot(R[i]) for i in range(R.shape[0])])

    betaE = np.vstack([beta(E[i]) for i in range(E.shape[0])])
    betaR = np.vstack([beta(R[i]) for i in range(R.shape[0])])
    COV = np.vstack(
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
    print COV

    #   perturbation: m = w+mbar ring elements
    p = SigmaP * np.vstack([ring_sample(s, 0, n) for i in range(m)])
