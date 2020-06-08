from prime import factorize, fast
from symbol_legendre import get_non_quad


def toneli_shanks(a, p):
    z = get_non_quad(p)
    q = p - 1
    s, q = factorize(q)
    c = fast(z, q, p)
    r = fast(a, (q + 1) // 2, p)
    t = fast(a, q, p)
    m = s
    while t % p != 1:
        i = 1
        while fast(t, 2 ** i, p) != 1:
            i += 1
        b = fast(c, 2 ** (m - i - 1), p)
        r = r * b % p
        t = t * b * b % p
        c = b * b % p
        m = i
    return r, p - r
