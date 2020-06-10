from prime import fast_sieve_initialize
from math import log2


# Returns true if number n is b-smooth. Uses division and modulus.
def brute_smooth(n, b):
    # primes = fast_sieve_initialize(b + 1)
    for p in primes:
        while n % p == 0:
            n //= p
    return n == 1 or n == -1


# Returns true if number n is b-smooth. Uses log2 and subtractions.
def is_b_smooth(n, b, primes, log_primes):
    log_counter = n.bit_length() - 1
    factors = [0] * len(primes)
    for p in range(0, len(primes)):
        while n % primes[p] == 0:
            log_counter -= log_primes[p]
            primes[p] **= 2
            factors[p] += 1
    if log_counter <= b.bit_length() - 1:
        return factors, True
    return factors, False


n1 = 216564934649977779183332104119550719684705171673087005634079598195092857334543
x1 = 1082881598137842789345394688  # 465365377579787109143940078957624295424
i1 = 1
b1 = 394840
primes = fast_sieve_initialize(b1)
print("init")
t1 = 0
while False:
    t1 = (i1 + x1) ** 2 - n1
    t2 = (x1 - i1) ** 2 - n1
    if brute_smooth(t1, b1):
        print("positive")
        break
    if brute_smooth(t2, b1):
        print("negative")
        break

    i1 += 1

print(t1)
