import numpy as np
import gaussian
import trapdoors
import hashlib

def hash_string(msg, m, q):
    h = hashlib.sha256(msg).digest()
    h_list = map(ord, h)
    x = 0
    for i in range(len(h_list)):
        x += 256**(i*q)*int(h_list[i])

    hash_array = []

    for j in range(m):
        hash_array.append(x % q)
        x /= q

    return np.array(hash_array)

class Signer:
    def __init__(self, n=8, m=150, sigma=4.7, q=2053):
        self.n, self.m, self.q = n,m,q
        self.A, self.R = trapdoors.gen_trap(n,q,m)
        self.B = trapdoors.gen_basis(n, q, m, self.A, self.R)
        self.sigma = sigma
        self.q = q

    def sign(self, msg):
        c = hash_string(msg, self.m, self.q)
        e = gaussian.gauss_samp(self.B, self.sigma, c, self.m, self.q)

        return c,e

class Verifier:
    def __init__(self, A, q=2053):
        self.A = A
        self.q = q

    def verify(self, msg, sig):
        (n, m) = self.A.shape
        h1 = hash_string(msg, m, self.q)

        h2 = np.mod( self.A*np.matrix(sig).T[0], self.q)

        print "h1: ", h1
        print "h2: ", h2

        return h1 == h2

if __name__ == "__main__":
    signer = Signer()

    msg = "hello world"
    c,e = signer.sign(msg)

    print "mean vector: ", c
    print "signature: ", e

    verifier = Verifier(signer.A)
    print "verified: ", verifier.verify(msg, e)
