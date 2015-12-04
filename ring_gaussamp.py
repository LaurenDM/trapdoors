import numpy as np
from ring_trapdoors import *
from gaussian import *
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
def ring_preimage_sample_G(k, u, s):
    x = np.zeros((k, len(u)))
    for i in range(len(u)):
        x[:,i] = preimage_sample_G_1D(k, u[i], s, len(u))
    return x

def ring_sample(sigma, mean, n):
    x = []
    for i in range(n):
        x.append(gauss_samp_1D(sigma, mean, n))
    return np.matrix(x)

def beta(x):
    return Rot(x)[0]

def preimage_sample_A(A, R, E, u, s, q, r):
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
    R1R2 = np.vstack([R1, R2, np.eye(R1.shape[1])])
    s = np.abs(np.linalg.svd(R1R2)[1][0]) * r * 32
    print s

    SigmaP = s**2*np.eye(COV.shape[0]) - COV
    rootSigma = np.linalg.cholesky(SigmaP - (r/2)**2*np.eye(SigmaP.shape[0]))

    #   perturbation: m = w+mbar ring elements
    p = rootSigma * np.matrix([np.random.normal(0, 1) for i in range(kplus2 * n)]).T
    randomizedRound = np.vectorize(lambda x: gauss_samp_1D(r/2, x, n))
    p = randomizedRound(p)

    p = np.array(p.reshape(kplus2, n))

    v = A_mult(q, A, p)
    preimage = ring_preimage_sample_G(k, u-v, r)
    return p + combine_sample(q, R, E, preimage)
