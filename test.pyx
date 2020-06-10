cimport cython
from libc.stdlib cimport malloc, free
from gmpy2 import xmpz, log2 as log_2


# Takes one number as argument and finds all primes up to that number with the sieve of eratosthenes.
@cython.boundscheck(False)
def fast_sieve_initialize(upper_bound):
    cdef char *sieve = <char *> malloc(upper_bound * sizeof(char))  # Initializing a dynamic array.
    if not sieve:  # In case of a memory error
        raise MemoryError()

    primes = [2]  # This is the final list that will hold the primes. Initializing with 2 in it as a special case.
    cdef long long int position = 3  # Start exploring from primes from 3.
    cdef long long int i  # iterator variable.
    cdef long long int p  # current prime variable.
    cdef long long int p_factor  # variable that is used to track all products of a prime.

    for i in range(0, upper_bound):
        sieve[i] = 0  # 0 means that the number is unexplored and possible prime
    i = 2
    try:
        while i < upper_bound:
            sieve[i] = 1  # Marking with 1 all numbers that can be divided by 2
            i += 2
        while position < upper_bound:  # Exploring all numbers below limit.
            if not sieve[position]:  # If we find a number that is 0, then it's a prime!
                p = position
                primes.append(p)
                p_factor = p
                while p_factor < upper_bound:  # Set with 1 all multiples of the prime
                    sieve[p_factor] = 1
                    p_factor += 2 * p  # Extra optimization! Set as 1 only the values
                                       # that are in form of p + 2*j for j>=0, because we already sieved multiples of 2.
            else:
                position += 2  # Just explore next odd number for prime.
    finally:
        free(sieve)  # Free up all that allocated space
    return primes

# These variables are used to find B-smooth numbers.
n = xmpz()
x = xmpz()
min_x = xmpz()
cdef long int z
cdef long int *a1
cdef long int *a2
cdef int *V
cdef long long int *p_list
cdef long int p_list_length
cdef int *log_primes

def initialize_sieve(py_n, py_x, py_z, py_a1, py_a2, py_list):
    global n  # The number to be factored.
    global x  # The ceil of root of n.
    global min_x  # Low-limit for searching b-smooth numbers.
    global z  # Range for searching b-smooth numbers. High-limit is defined as z + min_x.
    global a1 # toneli_shanks 1st solutions.
    global a2 # toneli_shanks 2nd solutions.
    global V  # Keeps track of possible primes.
    global p_list  # Prime number list
    global p_list_length  # Length of above list
    global log_primes  # Log2 of prime numbers
    #---Initialization---#
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
        log_primes[i] = int(round(log_2(py_list[i])))


# De-allocate all that space
def de_allocate_sieve():
    free(a1)
    free(a2)
    free(V)
    free(p_list)
    free(log_primes)


# The most interesting and annoying part... Finding B-smooth numbers!
def b_sieving():
    ret = []  # Return list of b-smooth numbers
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
    while not ret:  # Keep searching until there is a list to return.
        #Function to determine a low-bound of the logarithm of possible b-smooth numbers.
        smooth_limit = int(log_2((x + min_x) ** 2 - n)) - log_primes[-1] - 2

        # Sieving with 2 has been completely deleted as it didn't offer a lot...
        for j in range(1, p_list_length):
            # Sieve all values of a1[j] and its multiples.
            while a1[j] < z:
                V[a1[j]] += log_primes[j] # Add the approx. log2 of prime. (Additions are lot faster than divisions.)
                a1[j] += p_list[j] # Save the position of next multiple of the prime.
            a1[j] -= z # Subtract z to prepare for the next range.
            # Same procedure with a2[j] as above
            while a2[j] < z:
                V[a2[j]] += log_primes[j]
                a2[j] += p_list[j]
            a2[j] -= z
        for j in range(0, z):  # Check if any number is greater or equal than the smooth_limit.
            if V[j] >= smooth_limit:  # If it is, add it to list for further investigation.
                ret.append(int(x + j + min_x))
            V[j] = 0  # Set V for the next range of investigation.
        min_x += z  #Lower bound raises.
    return ret


# This function checks if a number is b_smooth and returns the list of powers of its primes factors.
# If not, it returns an empty list.
def factor(b_smooth_num, n):
    global p_list
    global p_list_length
    # check if B**2 - N is smooth.
    remainder = xmpz(b_smooth_num ** 2 - n)
    factor_list = [0] * p_list_length
    cdef long int j = 0
    # For each prime, trie to divide the remainder. Keep the successes in the factor_list for returning.
    while j < p_list_length:
        if remainder == 1:  # Check if remainder is 1, so return the list as all the factors were collected.
            return factor_list
        while remainder.__mod__(p_list[j]) == 0:
            remainder = remainder.__floordiv__(p_list[j])
            factor_list[j] += 1
        j += 1
    if remainder != 1:  # If remainder is not 1, then the number was not b-smooth. What a great sad ending!
        return []
    return factor_list