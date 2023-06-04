import contextlib
from millercheck import my_miller_rabin
from factorize import factorize
from typing import List, Tuple
from collections import Counter
from math import exp, log, e
from random import randint
from time import time
import math
import threading
import time



def calc_b_val(p: int) -> int:
    return round(exp(0.5 * pow(log(p, e) * log(log(p, e), e), 0.5)) * 3.38)


def generate_factor_base(b_val: int) -> List[int]:
    return [i for i in range(2, b_val) if my_miller_rabin(i)]


def generate_mtx(factor_base: List[int], p: int, a: int, overkill: int) -> Tuple[List[List[int]], List[int]]:
    fb_length = len(factor_base)
    n = p - 1
    mtx = []
    solutions_array = []

    while len(solutions_array) < fb_length + overkill:
        temp_vector, k = create_temp_vector(factor_base, n, a, p)
        if temp_vector is None:
            continue
        mtx.append(temp_vector)
        solutions_array.append(k)
    return mtx, solutions_array


def create_temp_vector(factor_base: List[int], n: int, a: int, p: int) -> Tuple[List[int] or None, int]:
    temp_vector = [0] * len(factor_base)
    k = randint(0, n - 1)
    a_k = (a**k) % p
    factors = factorize(a_k)
    if not valid_factors(factors, factor_base):
        return None, k

    return fill_temp_vector(factor_base, factors), k


def valid_factors(factors: List[int], factor_base: List[int]) -> bool:
    if not factors:
        return False
    return max(factors) <= factor_base[-1] and len(set(factors)) > 1


def fill_temp_vector(factor_base: List[int], factors: List[int]) -> List[int]:
    factor_counts = Counter(factors)
    return [factor_counts.get(b_val, 0) for b_val in factor_base]


def gaussian_elimination(mtx: List[List[int]], solutions_array: List[int], n: int):
    height, width = len(mtx), len(mtx[0])

    for j in range(height):
        for i in range(width):
            temp = mtx[j][i]
            if math.gcd(n, temp) != 1:
                continue
            rot = pow(temp, -1, n)
            adjust_rows(mtx, solutions_array, rot, j, i, n)
            break


def adjust_rows(mtx: List[List[int]], solutions_array: List[int], rot: int, j: int, i: int, n: int):
    height = len(mtx)
    
    mtx[j] = [(val * rot) % n for val in mtx[j]]
    solutions_array[j] = (solutions_array[j] * rot) % n

    for k in range(height):
        if k != j:
            coef = mtx[k][i]
            if coef != 0:
                adjust_row(mtx, solutions_array, coef, j, k, n)


def adjust_row(mtx: List[List[int]], solutions_array: List[int], coef: int, j: int, k: int, n: int):
    width = len(mtx[0])
    mtx[k] = [(mtx[k][wow] - mtx[j][wow] * coef) % n for wow in range(width)]
    solutions_array[k] = (solutions_array[k] - solutions_array[j] * coef) % n


def calc_result(mtx: List[List[int]], solutions_array: List[int], width: int) -> List[int]:
    result = [0] * width
    for i, row in enumerate(mtx):
        if 1 in row:
            with contextlib.suppress(ValueError):
                index = row.index(1)
                result[index] = solutions_array[i]
    return result


def find_res_log(factor_base: List[int], m: List[int], p: int, a: int, b: int, n: int) -> int:
    while True:
        l = randint(0, n - 1)
        b_a_l = b * pow(a, l, p) % p
        res_factors = factorize(b_a_l)
        if not valid_factors(res_factors, factor_base):
            continue
        res_factors = list(Counter(res_factors).items())
        res_log = calculate_res_log(res_factors, m, factor_base, l, n)
        if pow(a, res_log, p) == b:
            break
    return res_log


def calculate_res_log(res_factors: List[Tuple[int, int]], m: List[int], factor_base: List[int], l: int, n: int) -> int:
    res_log = -l
    for d in res_factors:
        index = factor_base.index(d[0])
        res_log += d[1] * m[index]
    res_log %= n
    return res_log


stop_threads = False
res_log = None

def generate_mtx_thread(factor_base: list, p: int, a: int, mtx: list, solutions_array: list, overkill: int):
    matr, sol_arr = generate_mtx(factor_base, p, a, overkill)  # Pass overkill parameter
    mtx.extend(matr)
    solutions_array.extend(sol_arr)
    
def stop_all_threads(res):
    global stop_threads
    global res_log
    stop_threads = True
    res_log = res
    stop_threads = False  # Reset the flag

def check_solutions_thread(a, b, p, factor_base, m, n):
    if not stop_threads:
        res = find_res_log(factor_base, m, p, a, b, n)
        if res is not None:
            stop_all_threads(res)





def main_calculus_threaded(a: int, b: int, p: int, threads=2) -> int or None:
    mtx, solutions_array = [], []
    b_val = calc_b_val(p)
    factor_base = generate_factor_base(b_val)
    fb_length = len(factor_base)
    overkill = int(fb_length * 0.1) + 7
    while not stop_threads:
        threads = []
        for _ in range(threads):
            t = threading.Thread(target=generate_mtx_thread, 
                                 args=(factor_base, p, a, mtx, solutions_array, overkill))
            threads.append(t)
            t.start()

        for t in threads:
            t.join()

        if len(solutions_array) < fb_length + overkill:
            continue
        
        n = p - 1
        gaussian_elimination(mtx, solutions_array, n)
        m = calc_result(mtx, solutions_array, fb_length + overkill)

        threads = []
        for _ in range(threads):
            t = threading.Thread(target=check_solutions_thread, args=(a, b, p, factor_base, m, n))
            threads.append(t)
            t.start()

        for t in threads:
            t.join()

        break

    return res_log

# thread_number = 2

# a = int(input("Введіть а: "))
# b = int(input("Введіть b: "))
# p = int(input("Введіть p: "))

# start_time = time.time()
# res_log = main_calculus_threaded(a, b, p)
# end_time = time.time()

# execution_time = end_time - start_time

# print(f'Час виконання: {execution_time} секунд')
# print(f"Результат: {res_log}")