
import numpy as np
import gaussian
import ring_gaussamp
import trapdoors
import ring_trapdoors
import hashlib
from numpy import log2, ceil


def hash_string(msg, n, q):
    h = hashlib.sha256(msg).digest()
    h_list = map(ord, h)
    x = 0
    for i in range(len(h_list)):
        x += (256**i)*int(h_list[i])

    hash_array = []

    for j in range(n):
        hash_array.append(x % q)
        x /= q

    return np.array(hash_array)

class Signer:
    def __init__(self, n=128, sigma=100000, k=19):
        self.n, self.k, = n,k
        self.q = 2**k
        self.r=9.4
        self.A, self.R, self.E = ring_trapdoors.gen_trap(n,self.q)
        self.rootSigma = ring_gaussamp.get_rootSigma(self.A,self.R,self.E,sigma, self.q,self.r)
        #self.B = trapdoors.gen_basis(n, q, m, self.A, self.R)
        self.sigma = sigma

    def privateKey(self):
        return (self.n, self.k, self.q, self.A, self.R, self.E, self.rootSigma, self.sigma)
    def publicKey(self):
        return (self.A, self.sigma, self.k)

    def sign(self, msg):
        u = hash_string(msg, self.n, self.q)

        #bin_dec = trapdoors.binary_decomp(np.matrix(u).T, int(ceil(log2(self.q))))
        #t = np.array(np.mod(self.R * bin_dec, self.q).T)[0]
        #v = gaussian.gauss_samp(self.B, self.sigma, -1*t, self.m, self.q)

        #sig = t+v

        sig = ring_gaussamp.preimage_sample_A(self.A,self.R,self.E,self.rootSigma, u,self.q,self.r)

        return sig

class Verifier:
    def __init__(self, A, sigma, k):
        self.A = A
        self.k = k
        self.q = 2**k
        self.sigma = sigma

    def verify(self, msg, sig):
        n = self.A.shape[1]

        # verify sig is small enough
        lengths = []
        for i in range(sig.shape[0]):
            lengths.append(np.linalg.norm(sig[i]))
        maxLength = np.max(lengths)

        # ensure that max norm vec of signature
        # is small enough, smaller than sigma*n*\sqrt(k)
        if  maxLength > self.sigma*n*np.sqrt(self.k):
            print("Signature too long: ", maxLength)
            print("Expecting signature smaller than ", self.sigma*n*np.sqrt(self.k))

            return False

        # generate the two hashes to compare
        h1 = hash_string(msg, n, self.q)
        h2 = np.mod( ring_trapdoors.A_mult(self.q,self.A,sig), self.q)

        # verify that h1 matches h2 everywhere
        return reduce(lambda x,y: x and y, (h1 == h2))
