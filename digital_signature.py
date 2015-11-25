import numpy as np
import gaussian
import trapdoors
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
    def __init__(self, n=128, m=2400, sigma=1000, q=2053):
        self.n, self.m, self.q = n,m,q
        self.A, self.R = trapdoors.gen_trap(n,q,m)
        self.B = trapdoors.gen_basis(n, q, m, self.A, self.R)
        self.sigma = sigma
        self.q = q

    def sign(self, msg):
        u = hash_string(msg, self.n, self.q)

        bin_dec = trapdoors.binary_decomp(np.matrix(u).T, int(ceil(log2(self.q))))
        t = np.array(np.mod(self.R * bin_dec, self.q).T)[0]
        v = gaussian.gauss_samp(self.B, self.sigma, -1*t, self.m, self.q)

        sig = t+v

        return sig

class Verifier:
    def __init__(self, A, q=2053):
        self.A = A
        self.q = q

    def verify(self, msg, sig):
        (n, m) = self.A.shape

        h1 = hash_string(msg, n, self.q)

        h2 = np.mod( np.array(self.A*np.matrix(sig).T).T[0], self.q)

        print "H1: ", h1
        print "H2: ", h2
        # verify that h1 matches h2 everywhere
        return reduce(lambda x,y: x and y, (h1 == h2))

if __name__ == "__main__":
    signer = Signer()

    msg = "hello world"
    e = signer.sign(msg)

    #print "signature: ", e

    verifier = Verifier(signer.A)

    print "verified: ", verifier.verify(msg, e)
