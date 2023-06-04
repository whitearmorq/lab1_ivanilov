import math


def floyd(n: int) -> int:

    def simple_trial_div(n: int) -> int:
        small_primes = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)
        return next((prime_number for prime_number in small_primes if n % prime_number == 0), n)
    
    def gcd(a: int, b: int) -> int:
        while b != 0:
            a, b = b, a % b
        return a
    
    tmp_factor = simple_trial_div(n)
    if tmp_factor != n:
        return tmp_factor
    
    for i in range(2, int(math.sqrt(n)) + 1):
        x_0 = i
        y_0 = i
        while True:
            x_0 = (pow(x_0, 2) + 1) % n
            y_0 = pow(((pow(y_0, 2) + 1) % n), 2) + 1
            y_0 = y_0 % n
            if x_0 == y_0:
                break
            d = gcd(abs(y_0 - x_0), n)
            if d != 1:
                return d
    
    return n


def factorize(p: int) -> list:
    result = []
    i = 2
    while pow(i, 2) <= p:
        if p % i:
            i += 1 if i == 2 else 2
        else:
            p //= i
            result.append(i)
    if p > 1:
        result.append(p)
    return result