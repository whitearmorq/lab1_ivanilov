import math

def f(x, n):
    return (x ** 2 + 1) % n

def pollards_rho(n):
    if n % 2 == 0:
        return 2
    x = 2; y = 2; d = 1
    while d == 1:
        x = f(x, n)
        y = f(f(y, n), n)
        d = math.gcd(abs(x - y), n)
    return d