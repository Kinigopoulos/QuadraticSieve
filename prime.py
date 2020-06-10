import random


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
