import random


# A fast algorithm which initializes eratosthenis's sieve.
def fast_sieve_initialize(upper_bound):
    sieve = [False] * upper_bound  # False means that the number is unexplored and possible prime
    primes = [2]
    position = 3
    for i in range(2, upper_bound, 2):
        sieve[i] = True  # Marking as True all numbers that can be divided by 2
    while position < upper_bound:
        if not sieve[position]:
            p = position
            primes.append(p)
            p_factor = p
            while p_factor < upper_bound:
                sieve[p_factor] = True
                p_factor += 2 * p
        else:
            position += 2
    return primes


# Checks if a number is divisible with numbers in [2, eratosthenis_max_number]
def eratosthenis_test(n):
    if n <= eratosthenis_max_number:
        return eratosthenis.__contains__(n)
    for p in eratosthenis:
        if n % p == 0:
            return False
    return True


# Returns the value of a^g mod n
def fast(a, g, n):
    x = a
    d = 1
    while g > 0:
        if g & 1:
            d = (d * x) % n
        g >>= 1
        x = (x * x) % n
    return d


# Returns true if it is a possible prime, using fermat's little theorem
def fermat_test(num):
    i = 2
    while i < 10:
        if fast(i, num - 1, num) != 1:
            return False
        i += 1
    return True


# Returns the biggest power of 2 which can factorize a number
def factorize(num):
    r = 0
    while num % 2 == 0:
        num //= 2
        r += 1
    return r, num


# Returns true if it is a prime with probability (1 - (1 / 4^k)) using Miller-Rabin test.
def miller_rabin_test(num, k):
    r, d = factorize(num)
    for i in range(0, k):
        a = random.randint(2, num - 2)
        x = fast(a, d, num)
        if x == 1 or x == num - 1:
            continue
        continue_loop = False
        for j in range(0, r - 1):
            x = fast(x, 2, num)
            if x == num - 1:
                continue_loop = True
                break
        if continue_loop:
            continue
        return True
    return False


# Sieve of Eratosthenis's Initialization
# Random Initialization
random.seed()
eratosthenis_max_number = 100000
eratosthenis = fast_sieve_initialize(eratosthenis_max_number)
