import numpy as np
from ring_trapdoors import *
def preimage_sample_G_1D(m, u, s):
    x = []
    for i in range(m):
        if u % 2 == 0:
            x_i = gauss_samp_1D_even()
        else:
            x_i = gauss_samp_1D_odd()
        x.append(x_i)
        u = (u-x_i) / 2
    return np.matrix(x)


#   m: sample count
#   u: target syndrome (ring element)
#   s: sigma
def ring_preimage_sample_G(m, u, s):
    x = []
    for i in range(len(u)):
        x.append(preimage_sample_G(m, u[i], s))
    return x

def ring_sample(sigma, mean, n):
    x = []
    for i in range(n):
        x.append(gauss_samp_1D(sigma, mean, n))
    return np.matrix(x)

def preimage_sample_A(A_0, R, u, s):
    print "A_0: " + str(np.shape(A_0))
    print A_0
    print "R:   " + str(np.shape(R))
    w = np.shape(R)[1]
    (n, mbar) = np.shape(A_0)
    m = w + mbar
    Sigma = s**2
    SigmaG = 2
    SigmaP = Sigma - np.concatenate([R, np.eye(w)]) * SigmaG * np.concatenate([R.T, np.eye(w)], axis=1)

    #   perturbation: m = w+mbar ring elements
    p = SigmaP * np.vstack([ring_sample(s, 0, n) for i in range(m)])

