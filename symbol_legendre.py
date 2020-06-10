from prime import fast


# Returns the (a|p) symbol.
def symbol(a, p):
    a %= p
    if a == 0:
        return 0
    result = fast(a, (p - 1) // 2, p)
    if result == 1:
        return result
    return -1


# Returns a number a so that (a|p)=-1
def get_non_quad(p):
    a = 2
    while a < p:
        if symbol(a, p) == -1:
            return a
        a += 1
    return 0
