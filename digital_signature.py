import numpy as np
import gaussian
import trapdoors
import hashlib
# class Signer:
#     def __init__(self, n, m, sigma, q):
#         self.A = np.zeros((n,m))
#         self.R = np.zeros((m,m))
#         self.sigma = sigma
#         self.q = q
#
#     def sign(self, m):
#         h = hashlib.sha256(m).hexdigest()


# outputs H(m) \in Z^n_q
def hash_message(m, n, q):
    h = hashlib.sha256(m).digest()

    x = 0
    for i in range(len(h)):
        x += (256**i)*h[i]


    hash_array = []
    for j in range(n):
        hash_array.append(x%q)
        x /= q

    return np.array(hash_array)

if __name__ == "__main__":
    n,q,m = 8,2053,150
    sigma = 4.7
    A,R = trapdoors.gen_trap(n,q,m)

    print "generated trapdoor: ", R
    B = trapdoors.gen_basis(n, q, m, A, R)

    message = "hello world"
    c = hash_message(message)

    print "mean vector: ", c

    e = gaussian.gauss_samp(B, sigma, c, n)

    print "signature: ", e
