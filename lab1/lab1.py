from collections import defaultdict
from math import gcd, isqrt, exp, log, sqrt
import random
import time


# Тест Міллера-Рабіна
def miller_rabin(p, k=10):
    if p <= 1:
        return False
    if p <= 3:
        return True
    if p % 2 == 0:
        return False

    d, s = p - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1

    def check_composite(a, d, p, s):
        x = pow(a, d, p)
        if x == 1 or x == p - 1:
            return False
        for _ in range(s - 1):
            x = pow(x, 2, p)
            if x == p - 1:
                return False
        return True

    for _ in range(k):
        a = random.randint(2, p - 2)
        if check_composite(a, d, p, s):
            return False
    return True


# Метод пробних ділень
def brute_force(n, limit=47):
    for i in range(2, min(n, limit) + 1):
        if n % i == 0:
            return i
    return None


# ρ-метод Полларда
def pollards_rho(n):
    if miller_rabin(n):
        return None

    def f(x, c, n):
        return (x * x + c) % n

    while True:
        x = y = 2
        c = random.randint(1, n - 1)
        d = 1

        while d == 1:
            x = f(x, c, n)
            y = f(f(y, c, n), c, n)
            d = gcd(abs(x - y), n)

        if d != 1 and d != n:
            return d


# Метод Брілхарта-Моррісона
# Побудова факторної бази
def build_factor_base(n, bound, k=10):
    base = [-1]
    for p in range(2, bound + 1):
        if miller_rabin(p, k) and pow(n, (p - 1) // 2, p) == 1:
            base.append(p)
    return base


# ланцюговий дріб sqrt n
def continued_fraction_sqrt(n, limit=200):
    u, v, a0 = 0, 1, isqrt(n)
    a = a0
    terms = [a0]
    for _ in range(limit):
        u = v * a - u
        v = (n - u * u) // v
        if v == 0:
            break
        a = (a0 + u) // v
        terms.append(a)
    return terms


# Перевірка B-гладкості
def is_b_smooth(n, base):
    factors = defaultdict(int)

    if n < 0:
        if -1 in base:
            factors[-1] = 1
        n = abs(n)

    for p in base:
        if p == -1:
            continue
        while n % p == 0:
            n //= p
            factors[p] += 1
    return n == 1, dict(factors)


def compute_L(n, a):
    ln_n = log(n)
    ln_ln_n = log(ln_n)
    return int(min(1000, exp(a * sqrt(ln_n * ln_ln_n))))


# метод Гауса
def gauss(matrix):
    m = len(matrix)
    n = len(matrix[0])
    A = [row[:] for row in matrix]
    cur_row = 0
    free_variables = []
    pivot_cols = []

    for col in range(n):
        pivot_row = -1
        for row in range(cur_row, m):
            if A[row][col] == 1:
                pivot_row = row
                break
        if pivot_row == -1:
            free_variables.append(col)
            continue

        A[cur_row], A[pivot_row] = A[pivot_row], A[cur_row]
        pivot_cols.append(col)

        for row in range(m):
            if row != cur_row and A[row][col] == 1:
                A[row] = [(A[row][i] ^ A[cur_row][i]) for i in range(n)]

        cur_row += 1

    basis = []
    for free_var in free_variables:
        solution = [0] * n
        solution[free_var] = 1

        for row in reversed(range(cur_row)):
            pivot_col = -1
            for col in range(n):
                if A[row][col] == 1:
                    pivot_col = col
                    break
            if pivot_col != -1 and pivot_col != free_var:
                s = sum(A[row][i] * solution[i] for i in range(pivot_col + 1, n)) % 2
                solution[pivot_col] = s

        basis.append(solution)

    return basis


def brillhart_morrison(n, k=10, max_a=6.0, step=0.25, time_limit=60):
    start_time = time.time()
    a = 1 / sqrt(2)

    while a <= max_a:
        if time.time() - start_time > time_limit:
            print(f"[{time.strftime('%H:%M:%S')}] Метод БМ не знайшов дільника за відведений час ({time_limit} с).")
            return None

        try:
            L = int(exp(a * sqrt(log(n) * log(log(n)))))
            base = build_factor_base(n, L, k)
        except Exception as e:
            print(f"[ERROR] Проблема з факторною базою: {e}")
            return None

        a_seq = continued_fraction_sqrt(n, 300)
        b1, b2 = 1, a_seq[1]
        bi_list = []
        exponent_vectors = []
        seen_vectors = {}

        for i in range(2, len(a_seq)):
            if time.time() - start_time > time_limit:
                print(f"[{time.strftime('%H:%M:%S')}] Метод БМ не знайшов дільника за відведений час.")
                return None

            bi = a_seq[i] * b2 + b1
            b1, b2 = b2, bi
            value = (bi * bi) % n

            smooth, factors = is_b_smooth(value, base)
            if not smooth:
                continue

            full_vector = [factors.get(p, 0) for p in base]
            mod2_vector = tuple(v % 2 for v in full_vector)

            if all(v == 0 for v in mod2_vector):  # костиль перший
                X, Y = bi % n, isqrt(value) % n
                for delta in [X - Y, X + Y]:
                    d = gcd(delta, n)
                    if 1 < d < n:
                        return d
                continue

            if mod2_vector in seen_vectors:  # костиль 2
                j = seen_vectors[mod2_vector]
                X = (bi * bi_list[j]) % n
                exp_sum = [full_vector[i] + exponent_vectors[j][i] for i in range(len(base))]
                Y_squared = 1
                for p, e in zip(base, exp_sum):
                    Y_squared *= pow(p, e // 2, n)
                    Y_squared %= n
                Y = isqrt(Y_squared) % n
                for delta in [X - Y, X + Y]:
                    d = gcd(delta, n)
                    if 1 < d < n:
                        return d
                continue

            seen_vectors[mod2_vector] = len(bi_list)
            bi_list.append(bi)
            exponent_vectors.append(full_vector)

        a += step

    print("[INFO] Жодного дільника не знайдено методом БМ")
    return None


def canonical_decomposition(n, k=10):
    factors = []
    start_time = time.time()
    pollard_used = False
    found_any = False  # Прапорець, чи хоч щось знайшли

    print(f"\n Початок розкладу числа {n} — {time.strftime('%H:%M:%S', time.localtime(start_time))}\n")

    while True:
        if miller_rabin(n, k):
            timestamp = time.strftime('%H:%M:%S', time.localtime())
            print(f"[{timestamp}] Просте число {n} додано до розкладу (Міллер-Рабін)")
            factors.append(n)
            found_any = True
            break

        divisor = brute_force(n)
        if divisor:
            timestamp = time.strftime('%H:%M:%S', time.localtime())
            print(f"[{timestamp}] Знайдено дільник {divisor} методом пробних ділень")
            factors.append(divisor)
            n //= divisor
            found_any = True
            continue

        if not pollard_used:
            divisor = pollards_rho(n)
            pollard_used = True
            if divisor:
                timestamp = time.strftime('%H:%M:%S', time.localtime())
                print(f"[{timestamp}] Знайдено дільник {divisor} ρ-методом Полларда")
                factors.append(divisor)
                n //= divisor
                found_any = True
                continue

        if miller_rabin(n, k):
            timestamp = time.strftime('%H:%M:%S', time.localtime())
            print(f"[{timestamp}] Просте число {n} додано до розкладу після Полларда")
            factors.append(n)
            found_any = True
            break

        divisor = brillhart_morrison(n, k)
        if divisor:
            timestamp = time.strftime('%H:%M:%S', time.localtime())
            print(f"[{timestamp}] Знайдено дільник {divisor} методом Брілхарта-Моррісона")
            factors.append(divisor)
            n //= divisor
            found_any = True
            continue
        else:
            break  # БМ не знайшов — виходимо

    end_time = time.time()
    print(f"\n Кінець роботи — {time.strftime('%H:%M:%S', time.localtime(end_time))}")

    if not found_any:
        print("\nНе вдалося знайти жодного дільника. Розклад завершено без результату.")
    else:
        print(f"\nКанонічний розклад числа: {sorted(factors)}")

    return sorted(factors)


# Тестування
number = 59
canonical_decomposition(number)



