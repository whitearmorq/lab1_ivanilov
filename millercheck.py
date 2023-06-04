from random import randrange

def calculate_d_and_s(n):
    s = 0
    d = n - 1

    while d % 2 == 0:
        d //= 2
        s += 1

    return s, d

def perform_check(a, s, d, n):
    x = pow(a, d, n)
    if x in [1, n - 1]:
        return True
    for _ in range(s - 1):
        x = pow(x, 2, n)
        if x == n - 1:
            return True
    return False

def my_miller_rabin(n, k=5):  
    if n in [2, 3]:
        return True

    if n <= 1 or n % 2 == 0:
        return False

    s, d = calculate_d_and_s(n)

    for _ in range(k):
        a = randrange(2, n - 1)
        if not perform_check(a, s, d, n):
            return False
    return True
