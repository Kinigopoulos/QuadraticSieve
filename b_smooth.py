from test import fast_sieve_initialize


# Returns true if number n is b-smooth. Uses division and modulus.
def brute_smooth(n, b):
    primes = fast_sieve_initialize(b + 1)
    for p in primes:
        while n % p == 0:
            n //= p
    return n == 1 or n == -1
