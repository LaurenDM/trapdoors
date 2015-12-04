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

    def sign(self, msg):
        u = hash_string(msg, self.n, self.q)

        #bin_dec = trapdoors.binary_decomp(np.matrix(u).T, int(ceil(log2(self.q))))
        #t = np.array(np.mod(self.R * bin_dec, self.q).T)[0]
        #v = gaussian.gauss_samp(self.B, self.sigma, -1*t, self.m, self.q)

        #sig = t+v

        sig = ring_gaussamp.preimage_sample_A(self.A,self.R,self.E,self.rootSigma, u,self.q,self.r)

        return sig

class Verifier:
    def __init__(self, A, k=19):
        self.A = A
        self.q = 2**k

    def verify(self, msg, sig):
        n = self.A.shape[1]

        h1 = hash_string(msg, n, self.q)

        h2 = np.mod( ring_trapdoors.A_mult(self.q,self.A,sig), self.q)

        print "H1: ", h1
        print "H2: ", h2
        # verify that h1 matches h2 everywhere
        return reduce(lambda x,y: x and y, (h1 == h2))

if __name__ == "__main__":
    signer = Signer(n=256)

    msg = "hello world"
    e = signer.sign(msg)

    print "\nsignature: ", e

    verifier = Verifier(signer.A)

    print "\nverified: ", verifier.verify(msg, e)
