from sage.all import *


def solve(a, b, mod):
    try:
        A = matrix(GF(mod), a).transpose()
        b = matrix(GF(mod), b)
        return A.right_kernel(b).list()
    except ValueError:
        return -1


def QS_Sieve(num):
    x = qsieve(num)
    print(x)
