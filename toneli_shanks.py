from prime import factorize, fast
from symbol_legendre import get_non_quad


# Given numbers a and p, this algorithm returns a number x so that (x ^ 2) mop p = a
def toneli_shanks(a, p):
# find z such as Legendre symbol (z/p) == -1
    z = get_non_quad(p)
    q = p - 1
    s, q = factorize(q) #(p-1)=2^n*k
    c = fast(z, q, p)
    r = fast(a, (q + 1) // 2, p)
    t = fast(a, q, p)
    m = s
    while t % p != 1:
        i = 1 #Find the least i such as t^2^imodp=1modp
        while fast(t, 2 ** i, p) != 1: #if i!=0 then increase until c^2^(m-i-1) = 1modp
            i += 1
        b = fast(c, 2 ** (m - i - 1), p)
        r = r * b % p
        t = t * b * b % p
        c = b * b % p
        m = i
    return r, p - r #solutions
