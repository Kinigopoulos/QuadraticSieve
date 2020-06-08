from prime import fast_sieve_initialize, fermat_test, miller_rabin_test, eratosthenis_test
from b_smooth import is_b_smooth, find_b_smooth_numbers, find_b_smooth_numbers2
import symbol_legendre
from modulo_equation import solve
from toneli_shanks import toneli_shanks
from math import exp, log, ceil, gcd, isqrt


def factor_base(n):
    if n > 10 ** 13:
        log_n = log(n)
        e_x = (log_n * log(log_n)) ** 0.5
        e = exp(e_x)
        return int(ceil(e ** 0.54))
    return int(n ** 0.5)


def QuadraticSieve(n):
    # Initialization: B-smooth number selection and prime list until B.
    B = factor_base(n)
    print("B:", B)
    # return 1
    prime_list = fast_sieve_initialize(B + 1)

    # Keep primes that are quadratic residue.
    p_list = [2]
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
            while remainder % p_list[j] == 0:
                remainder //= p_list[j]
                factor_list[j] += 1
        return factor_list

    z = 50000000
    pos = 0
    eq = []
    smooth_numbers = []
    final_factor = 1
    tries = 0
    sol = []
    while final_factor == 1 or final_factor == n:
        tries += 1
        print(len(sol))
        print("Finding smooth numbers.", tries, "tries. This procedure will take the longest...", end=" ")
        while True:
            S, pos, a1, a2 = find_b_smooth_numbers(n, x, z, p_list, a1, a2, pos)

            print("S:", S)

            for s in S:
                smooth_numbers.append(s)
                eq.append(factor(s))

                sol = solve(eq, [0] * S.__len__(), 2)

                if len(sol) > tries:
                    print("Solutions:", sol)
                    print(smooth_numbers)
                    break
            if len(sol) > tries:
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
        print(n, non_factorised_numbers)
        non_factorised_numbers.pop(0)
        if n == 1:
            continue
        elif n < eratosthenis_upper_bound:
            final_factors.append(n)
            continue
        elif fermat_test(n) and miller_rabin_test(n, 10):
            print(eratosthenis_test(n))
            final_factors.append(n)
            continue

        if not eratosthenis_test(n):
            print("eratos", n)
        elif not fermat_test(n):
            print("fermat", n)
        elif not miller_rabin_test(n, 10):
            print("miller", n)

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


eratosthenis_upper_bound = 150
small_factors = fast_sieve_initialize(eratosthenis_upper_bound)
print("Initializing sieve: ")
# QuadraticSieve(15347)
# print_factors(53987325434578734538546895332)

# res = QuadraticSieve((2 ** 277 - 1) // 1121297)

# print(2 ** 277 - 1)
# print(factor_base(2 ** 277 - 1))

# QuadraticSieve(5398739)
print_factors(2 ** 277 - 1)
# print_factors(2 ** 277 - 1)
