cimport cython
from libc.stdlib cimport malloc, free
from gmpy2 import xmpz, log2


@cython.boundscheck(False)
def fast_sieve_initialize(upper_bound):
    cdef char *sieve = <char *> malloc(upper_bound * sizeof(char)) # False means that the number is unexplored and possible prime
    if not sieve:
        raise MemoryError()

    primes = [2]
    cdef long long int position = 3
    cdef long long int i
    cdef long long int p
    cdef long long int p_factor

    for i in range(0, upper_bound):
        sieve[i] = 0
    i = 2
    try:
        while i < upper_bound:
            sieve[i] = 1  # Marking as True all numbers that can be divided by 2
            i += 2
        while position < upper_bound:
            if not sieve[position]:
                p = position
                primes.append(p)
                p_factor = p
                while p_factor < upper_bound:
                    sieve[p_factor] = 1
                    p_factor += 2 * p
            else:
                position += 2
    finally:
        free(sieve)
    return primes

n = xmpz()
x = xmpz()
min_x = xmpz()
cdef long int z
cdef long int *a1
cdef long int *a2
cdef int *V
cdef long long int *p_list
cdef int p_list_length
cdef int *log_primes

def initialize_sieve(py_n, py_x, py_z, py_a1, py_a2, py_V, py_list):
    global n
    global x
    global min_x
    global z
    global a1
    global a2
    global V
    global p_list
    global p_list_length
    global log_primes
    n = xmpz(py_n)
    x = xmpz(py_x)
    min_x = xmpz(0)
    z = py_z
    cdef long int length = len(py_a1)
    a1 = <long int *> malloc(length * sizeof(long int))
    a2 = <long int *> malloc(length * sizeof(long int))
    V = <int *> malloc(z * sizeof(int))

    p_list_length = len(py_list)
    p_list = <long long int *> malloc(p_list_length * sizeof(long long int))
    log_primes = <int *> malloc(p_list_length * sizeof(int))
    cdef long int i
    for i in range(0, z):
        V[i] = 0
    for i in range(0, length):
        a1[i] = py_a1[i]
        a2[i] = py_a2[i]
    i = 0
    for i in range(0, p_list_length):
        p_list[i] = py_list[i]
        log_primes[i] = int(round(log2(py_list[i])))


def de_allocate_sieve():
    free(a1)
    free(a2)
    free(V)
    free(p_list)
    free(log_primes)


def b_sieving():
    ret = []
    global n
    global x
    global min_x
    global z
    global a1
    global a2
    global V
    global p_list
    global p_list_length
    global log_primes
    cdef long int j
    cdef int smooth_limit
    while not ret:
        smooth_limit = int(log2((x + min_x) ** 2 - n)) - log_primes[-1] - 2

        for j in range(1, p_list_length):
            while a1[j] < z:
                V[a1[j]] += log_primes[j]
                a1[j] += p_list[j]
            a1[j] -= z
            while a2[j] < z:
                V[a2[j]] += log_primes[j]
                a2[j] += p_list[j]
            a2[j] -= z
        for j in range(0, z):
            # print(V[j], " <", smooth_limit)
            if V[j] >= smooth_limit:
                ret.append(int(x + j + min_x))
            V[j] = 0
        min_x += z
        # print(n, x, min_x, z)
    return ret


def factor(b_smooth_num, n, p_list):
    cdef long int length = len(p_list)
    remainder = xmpz(b_smooth_num ** 2 - n)
    factor_list = [0] * length
    cdef long int j = 0
    while j < length:
        if remainder == 1:
            break
        while remainder.__mod__(p_list[j]) == 0:
            remainder = remainder.__floordiv__(p_list[j])
            factor_list[j] += 1
        j += 1
    if remainder != 1:
        return []
    return factor_list