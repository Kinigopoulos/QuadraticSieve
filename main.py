from prime import fermat_test, miller_rabin_test
import symbol_legendre
from modulo_equation import solve
from toneli_shanks import toneli_shanks
from math import exp, log, ceil, gcd, isqrt
from timeit import default_timer as timer
from test import b_sieving, factor, initialize_sieve, fast_sieve_initialize, de_allocate_sieve


# Computes and returns B, given N.
def factor_base(n):
    if n > 10 ** 13:
        log_n = log(n)
        e_x = (log_n * log(log_n)) ** 0.5
        e = exp(e_x)
        return int(ceil(e) ** 0.5)
    return int(n ** 0.35)


def QuadraticSieve(n):
    # Initialization: B-smooth number selection and prime list until B.
    B = factor_base(n)
    print("B:", B)

    prime_list = fast_sieve_initialize(B + 1)  # Get a prime list with prime numbers up to B.

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
    # Find (a(i) ^ 2) mod 2 = n with one if. Also, 2 has one solutions so it's stored on a1.
    if (x ** 2 - n) % 2 != 0:
        a1[0] = 1
    print("Collecting appropriate primes from the sieve...", end="\t")
    for p in p_list[1:]:  # Execute toneli shanks with all primes, excluding 2.
        x1, x2 = toneli_shanks(n, p)
        if x1 >= 0:  # In case of a negative result, it means that the algorithm found no solution.
            a1.append((x1 - x) % p)
            a2.append((x2 - x) % p)
        else:
            p_list.remove(p)  # Just in case there are no solutions, remove the prime from the list.
    print("DONE")

    # Sieving part. Starting from the ceil of root of N to find B-Smooth numbers.
    z = p_list[-1] * 3  # Range for sieving B-smooth numbers.
    equations = []  # Stores all relations exported by factor(). More on that below.
    smooth_numbers = []  # Saves all B-smooth numbers.
    final_factor = 1  # The saint-number we are going to return, once finished.
    tries = 0  # Trial and retry. In case of failure, we'll look for more b-smooth numbers and solutions.
    solutions = []  # Solutions exported from modulo_equation().

    # Initialization of sieve with Cython.
    initialize_sieve(n, x, z, a1, a2, p_list)

    # While our final_factor is not an actual factor, keep TRYING!!!
    while final_factor == 1 or final_factor == n:
        tries += 1
        print(len(solutions))
        print("Finding smooth numbers.", tries, "tries. This procedure will take the longest...", end=" ")
        while True:
            S = b_sieving()  # Get the possible B-smooth numbers.
            print("Found these smooth numbers:", S, "Completion:", len(smooth_numbers), "/", len(p_list))

            for s in S:
                factor_matrix = factor(s, n)  # Factor each possible B-smooth number and make sure it's smooth.
                if factor_matrix:  # For each B-smooth number found
                    smooth_numbers.append(s)  # KEEP IT, they're precious!
                    equations.append(factor_matrix)  # Append the list exported from factors, they will be needed for
                    # the next step.

                    # If we're close to finding relations, start searching.
                    if len(smooth_numbers) >= len(p_list) * 0.95:
                        solutions = solve(equations, 2)
                        print("Checking for solutions from the modulus equation...")
                        # solutions[0] will always give an array filled with zeros and we don't want that.
                        if solutions and len(solutions) > tries:
                            print("Solutions:", solutions)
                            print(smooth_numbers)
                            break
            if solutions and len(solutions) > tries:
                break

        print("DONE")
        num_sol = 1
        while num_sol < len(solutions):
            sol2 = []
            for i in solutions[num_sol]:
                sol2.append(int(i))
            x1 = 1
            x2 = 1
            for i in range(0, sol2.__len__()):
                if sol2[i] > 0:
                    x1 *= (smooth_numbers[i] ** 2 - n) ** sol2[i]  # Define x1.
                    x2 *= (smooth_numbers[i] ** sol2[i]) ** 2      # Define x2.
            # Get their square roots.
            x1 = isqrt(x1)
            x2 = isqrt(x2)

            # Find their their gcd with n, by adding them or subtracting them.
            # If the value returned is not 1 or n, congratulate yourself, you just found a factor!
            final_factor = gcd(x2 + x1, n)
            print("___SOLUTION___: ", final_factor)
            if final_factor != 1 and final_factor != n:
                break
            final_factor = gcd(x2 - x1, n)
            print("___SOLUTION___: ", final_factor)
            if final_factor != 1 and final_factor != n:
                break
            num_sol += 1
    de_allocate_sieve()  # Deallocate all that memory from Cython's part.
    print(final_factor)
    return final_factor  # Return the prize!


# This function check if the number n has a factor in the eratosthenes's sieve and returns it.
# Otherwise, it invokes Quadratic Sieve.
def get_factor(n):
    for p in small_factors:
        if n % p == 0:
            return p
    return QuadraticSieve(n)


# Collects all factors of the number. Stops until all numbers in the list are primes.
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
            final_factors.append(n)
            continue

        factor1 = get_factor(n)

        factor2 = n // factor1
        non_factorised_numbers.append(factor1)
        non_factorised_numbers.append(factor2)

    return final_factors


# Prints beautifully all the factors of the number n.
def print_factors(n):
    non_factorised_numbers = [n]
    final_factors = factorise(non_factorised_numbers)
    print()
    print(n, "'s factors are:")
    for f in final_factors[:len(final_factors) - 1]:
        if f != 1:
            print(f, end=" x ")
    if final_factors[-1] != 1:
        print(final_factors[-1])


# Eratosthenes's sieve limit. Change it to taste.
eratosthenis_upper_bound = 100000000
small_factors = fast_sieve_initialize(eratosthenis_upper_bound)
print("Initializing sieve: ")

start = timer()
print_factors(2 ** 223 - 1)  # Insert your number to factorise here.
end = timer()

print("Completed in", end - start, "seconds")
