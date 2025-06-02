import sympy
import random
import time
import multiprocessing
from math import gcd


def factor_base(B):
    return list(sympy.primerange(2, B + 1))


def find_one_relation(p, g, FB):
    p_minus_1 = p - 1
    while True:
        e = random.randrange(0, p_minus_1)
        val = pow(g, e, p)
        if val == 1:
            continue
        fac = sympy.factorint(val)
        if all(prime in FB for prime in fac):
            exp_vector = [fac.get(prime, 0) for prime in FB]
            return e, exp_vector


def solve_linear_system(A, b, mod):
    n = len(A)
    m = len(A[0])
    M = [row[:] + [b_i] for row, b_i in zip(A, b)]
    row = 0
    for col in range(m):
        if row >= n:
            break
        pivot = None
        for r in range(row, n):
            if M[r][col] % mod != 0 and gcd(M[r][col], mod) == 1:
                pivot = r
                break
        if pivot is None:
            continue
        M[row], M[pivot] = M[pivot], M[row]
        inv = sympy.mod_inverse(M[row][col], mod)
        for c in range(col, m + 1):
            M[row][c] = (M[row][c] * inv) % mod
        for r in range(n):
            if r != row and M[r][col] != 0:
                factor = M[r][col]
                for c in range(col, m + 1):
                    M[r][c] = (M[r][c] - factor * M[row][c]) % mod
        row += 1
    x = [0] * m
    for r in range(n):
        for c in range(m):
            if M[r][c] == 1:
                x[c] = M[r][m] % mod
                break
    return x


def index_calculus(p, g, h, parallel=False, B=None):
    if B is None:
        c = 3.38
        logn = sympy.log(p - 1)
        loglogn = sympy.log(logn)
        B = int(c * sympy.exp(0.5 * (logn * loglogn) ** 0.5))

    FB = factor_base(B)
    r = len(FB)
    relations = []
    eq_matrix = []
    p_minus_1 = p - 1

    def collect_relations_serial():
        while len(eq_matrix) < r:
            e, exp_vec = find_one_relation(p, g, FB)
            relations.append((e, exp_vec))
            eq_matrix.append(exp_vec)

    def collect_relations_parallel():
        with multiprocessing.Pool() as pool:
            while len(eq_matrix) < r:
                to_do = min(multiprocessing.cpu_count(), r - len(eq_matrix))
                tasks = [pool.apply_async(find_one_relation, (p, g, FB)) for _ in range(to_do)]
                for task in tasks:
                    e, exp_vec = task.get()
                    if exp_vec not in eq_matrix:
                        relations.append((e, exp_vec))
                        eq_matrix.append(exp_vec)

    if parallel:
        collect_relations_parallel()
    else:
        collect_relations_serial()

    A = [vec for (_, vec) in relations]
    b = [e for (e, _) in relations]
    logs = solve_linear_system(A, b, p_minus_1)

    while True:
        k = random.randrange(0, p_minus_1)
        val = (pow(g, k, p) * h) % p
        if val == 1:
            continue
        fac = sympy.factorint(val)
        if all(prime in FB for prime in fac):
            f_vec = [fac.get(prime, 0) for prime in FB]
            x = sum(f_vec[i] * logs[i] for i in range(r)) - k
            return x % p_minus_1


def test_index_calculus():
    bit_sizes = [38, 42, 8, 20, 48]
    for bits in bit_sizes:
        p = sympy.nextprime(2 ** bits)
        g = sympy.primitive_root(p)
        x_true = random.randrange(1, p - 1)
        h = pow(g, x_true, p)

        times_serial = []
        for _ in range(5):
            t0 = time.time()
            x_found = index_calculus(p, g, h, parallel=False)
            dt = time.time() - t0
            times_serial.append(dt)
            if dt > 300:
                break
        avg_serial = sum(times_serial) / len(times_serial)

        times_parallel = []
        for _ in range(5):
            t0 = time.time()
            x_found = index_calculus(p, g, h, parallel=True)
            dt = time.time() - t0
            times_parallel.append(dt)
            if dt > 300:
                break
        avg_parallel = sum(times_parallel) / len(times_parallel)

        print(f"p (bits={bits}): avg_serial={avg_serial:.5f}s, avg_parallel={avg_parallel:.5f}s")


if __name__ == "__main__":
    test_index_calculus()
