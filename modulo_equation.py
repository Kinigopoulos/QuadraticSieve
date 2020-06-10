from sage.all import *


# Given a 2d-array A, give 1d-arrays X so that A*X=0 (mod 2)
def solve(a, mod):
    try:
        ring = IntegerModRing(mod)
        A = matrix(ring, a).transpose()

        return A.right_kernel().list()
    except ValueError:
        print("Value Error from Sage")
        return []


# This is the QS_Sieve implementation from SageMath.
def QS_Sieve(num):
    x = qsieve(num)
    print(x)
