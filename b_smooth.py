from prime import fast_sieve_initialize
from math import log2
from timeit import default_timer as timer


def is_b_smooth(n, b):
    primes = fast_sieve_initialize(b + 1)
    log_counter = log2(n)
    for p in primes:
        if n % p == 0:
            log_counter -= log2(p)
    if log_counter <= log2(b):
        return True
    return False


def find_b_smooth_numbers(n, b, z, primes, a1, a2, pos):
    initial_position = pos
    V = [float(0)] * z
    Y = [0] * z
    smooth_numbers = []

    while len(smooth_numbers) < 1:

        def log_2(x):
            return x.bit_length() - 1

        print(initial_position)
        for i in range(0, z):
            num = ((b + i + initial_position) ** 2) - n
            Y[i] = num
            V[i] = log_2(num)

        if Y[0] % 2 == 0:
            position = 0
        else:
            position = 1
        while position < len(V):
            V[position] -= 1
            position += 2

        for i in range(1, len(primes)):
            while a1[i] < len(V) + initial_position:
                V[a1[i] - initial_position] -= log_2(primes[i])
                a1[i] += primes[i]
            while a2[i] < len(V) + initial_position:
                V[a2[i] - initial_position] -= log_2(primes[i])
                a2[i] += primes[i]

        for i in range(0, len(V)):
            if V[i] <= log_2(primes[-1]):
                smooth_numbers.append(b + i + initial_position)

        initial_position += z
    return smooth_numbers, initial_position, a1, a2


def find_b_smooth_numbers2(n, b, z, primes, a1, a2, pos):
    initial_position = pos

    def sieve(x, divisor):
        position = x
        while position < len(V):
            while V[position] % divisor == 0:
                V[position] //= divisor
            position += divisor

    smooth_numbers = []
    while len(smooth_numbers) < len(primes):
        print(initial_position)
        Y = []
        V = []
        for i in range(0, z):
            num = ((b + i + initial_position) ** 2) - n
            Y.append(num)
            V.append(num)

        if Y[0] % 2 == 0:
            sieve(0, 2)
        else:
            sieve(1, 2)

        for i in range(1, len(primes)):
            sieve(a1[i], primes[i])
            sieve(a2[i], primes[i])

        for i in range(0, len(V)):
            if V[i] <= 1:
                smooth_numbers.append(b + i)

        initial_position += z
    return smooth_numbers, pos, a1, a2
