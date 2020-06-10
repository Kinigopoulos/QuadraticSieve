from sage.all import *


def solve(a, mod):
    try:
        ring = IntegerModRing(mod)
        A = matrix(ring, a).transpose()

        return A.right_kernel().list()
    except ValueError:
        print("Value Error from Sage")
        return []


def QS_Sieve(num):
    x = qsieve(num)
    print(x)
