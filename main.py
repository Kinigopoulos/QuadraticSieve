from prime import fast_sieve_initialize, fermat_test, miller_rabin_test, eratosthenis_test
from b_smooth import is_b_smooth, brute_smooth
import symbol_legendre
from modulo_equation import solve
from toneli_shanks import toneli_shanks
from math import exp, log, ceil, gcd, isqrt
from timeit import default_timer as timer
from test import b_sieving as c_sieve, factor as c_factor, initialize_sieve, fast_sieve_initialize as eratosthenis


def factor_base(n):
    if n > 10 ** 13:
        log_n = log(n)
        e_x = (log_n * log(log_n)) ** 0.5
        e = exp(e_x)
        return int(ceil(e) ** 0.525)
    return int(n ** 0.35)


def QuadraticSieve(n):
    # Initialization: B-smooth number selection and prime list until B.
    B = factor_base(n)
    print("B:", B)

    prime_list = eratosthenis(B + 1)  # Get a prime list with prime numbers up to B.


    # Keep primes that are quadratic residue.
    p_list = [2]  # 2 will always be in the list.
    for p in prime_list[1:]:
        if symbol_legendre.symbol(n, p) == 1:
            p_list.append(p)
    print("p:", p_list)
    print(len(p_list))

    # Find roots a(i) for each prime p => (a(i) ^ 2) mod p = n
    x = int(ceil(n ** 0.5))
    print("x:", x)
    a1 = [0]
    a2 = [0]
    print("Collecting appropriate primes from the sieve...", end="\t")
    for p in p_list[1:]:  # prime number 2 is a special occasion and it will be reviewed later
        x1, x2 = toneli_shanks(n, p)
        if x1 >= 0:  # In case of a negative result, it means that the algorithm found no solution.
            a1.append((x1 - x) % p)
            a2.append((x2 - x) % p)
        else:
            p_list.remove(p)
    print("DONE")

    # Sieving part. Starting from the ceil of root of N to find B-Smooth numbers.
    # This function will help later... It returns a list with the factors of a b-smooth number.
    def factor(b_smooth_num):
        remainder = b_smooth_num ** 2 - n
        factor_list = [0] * len(p_list)
        for j in range(0, len(p_list)):
            if remainder == 1:
                break
            while remainder % p_list[j] == 0:
                remainder //= p_list[j]
                factor_list[j] += 1
        if remainder != 1:
            return []
        return factor_list

    z = p_list[-1] * 3
    pos = 0
    eq = []
    smooth_numbers = []
    final_factor = 1
    tries = 0
    sol = []

    V = [0] * z
    if (x ** 2 - n) % 2 != 0:
        a1[0] = 1

    log_primes = p_list.copy()
    for i in range(0, len(log_primes)):
        log_primes[i] = log_primes[i].bit_length() - 1

    initialize_sieve(n, x, z, a1, a2, V, p_list)
    while final_factor == 1 or final_factor == n:
        def b_sieving(min_x):
            ret = []
            while len(ret) < 1:
                smooth_limit = ((x + min_x) ** 2 - n).bit_length() - log_primes[-1] - 2

                for j in range(1, len(p_list)):
                    while a1[j] < len(V) + min_x:
                        V[a1[j] - min_x] += log_primes[j]
                        a1[j] += p_list[j]
                    while a2[j] < len(V) + min_x:
                        V[a2[j] - min_x] += log_primes[j]
                        a2[j] += p_list[j]

                for j in range(0, len(V)):
                    # print(V[j], "<", smooth_limit, len(smooth_numbers))
                    if V[j] >= smooth_limit:
                        ret.append(x + j + min_x)
                    V[j] = 0
                min_x += z
            return ret, min_x

        tries += 1
        print(len(sol))
        print("Finding smooth numbers.", tries, "tries. This procedure will take the longest...", end=" ")
        while True:
            # S, pos = b_sieving(pos)  # find_b_smooth_numbers(n, x, z, p_list, a1, a2, pos)
            S = c_sieve()
            # break
            print("Found these smooth numbers:", S, "Completion:", len(smooth_numbers), "/", len(p_list))

            for s in S:
                # factor_matrix, is_smooth = b_smooth(s ** 2 - n, B, p_list)
                #if not is_smooth:
                #    continue
                # print("SMOOTHNESS:", brute_smooth(s**2-n, B + 1), "///", s)
                factor_matrix = c_factor(s, n, p_list)
                if factor_matrix:
                    smooth_numbers.append(s)
                    eq.append(factor_matrix)
                    if len(smooth_numbers) >= len(p_list):
                        sol = solve(eq, [0] * len(smooth_numbers), 2)
                        if sol and len(sol) > tries:
                            print("Solutions:", sol)
                            print(smooth_numbers)
                            break
            if sol and len(sol) > tries:
                break

        print("DONE")
        num_sol = 1
        while num_sol < len(sol):
            sol2 = []
            for i in sol[num_sol]:
                sol2.append(int(i))
            x1 = 1
            x2 = 1
            print(sol2)
            for i in range(0, sol2.__len__()):
                if sol2[i] > 0:
                    x1 *= (smooth_numbers[i] ** 2 - n) ** sol2[i]
                    x2 *= (smooth_numbers[i] ** sol2[i]) ** 2

            x1 = isqrt(x1)
            x2 = isqrt(x2)

            final_factor = gcd(x2 - x1, n)
            print("___SOLUTION___: ", final_factor)
            final_factor = gcd(x2 + x1, n)
            print("___SOLUTION___: ", final_factor)
            if final_factor != 1 and final_factor != n:
                break
            num_sol += 1
    print(final_factor)
    print(pos)
    return final_factor


def get_factor(n):
    for p in small_factors:
        if n % p == 0:
            return p
    return QuadraticSieve(n)


def factorise(non_factorised_numbers):
    final_factors = []
    while non_factorised_numbers:
        n = non_factorised_numbers[0]
        print("Factoring", n, "with these numbers left on the list", non_factorised_numbers)
        non_factorised_numbers.pop(0)
        if n == 1:
            continue
        elif eratosthenis_upper_bound > n == get_factor(n):
            final_factors.append(n)
            continue
        elif fermat_test(n) and miller_rabin_test(n, 10):
            print(eratosthenis_test(n))
            final_factors.append(n)
            continue

        factor1 = get_factor(n)

        factor2 = n // factor1
        non_factorised_numbers.append(factor1)
        non_factorised_numbers.append(factor2)

    return final_factors


def print_factors(n):
    non_factorised_numbers = [n]
    final_factors = factorise(non_factorised_numbers)
    print(n, "'s factors are:")
    for f in final_factors[:len(final_factors) - 1]:
        if f != 1:
            print(f, end=" x ")
    if final_factors[-1] != 1:
        print(final_factors[-1])


eratosthenis_upper_bound = 100000000
small_factors = eratosthenis(eratosthenis_upper_bound)
print("Initializing sieve: ")

start = timer()
print_factors(2 ** 277 - 1)
end = timer()

print("Completed in", end - start, "seconds")
