import time
from collections import defaultdict


def euclid_gcd(a, b):
    while b:
        a, b = b, a % b
    return a
 
def create_cyclic_group(p: int):
    return [i for i in range(1, p) if euclid_gcd(i, p) == 1]


def extended_euclid_gcd(a, b):
    if a == 0:
        return b, 0, 1
    g, y, x = extended_euclid_gcd(b % a, a)
    return g, x - (b // a) * y, y

def inverse_mod(a, mod):
    g, x, _ = extended_euclid_gcd(a, mod)
    if g != 1:
        raise ValueError("Оберненого елемента не існує")
    return x % mod


def dissolve(n):
    divs = defaultdict(int)
    for d in range(2, int(n/2) + 1):
        while n % d == 0:
            divs[d] += 1
            n //= d
    if n > 1:
        divs[n] += 1
    return dict(divs)


def power_mod(base, exponent, mod):
    result = 1
    base %= mod
    while exponent:
        if exponent & 1:
            result = (result * base) % mod
        base = (base * base) % mod
        exponent >>= 1
    return result


def chinese_remainder_theorem(pairs):
    result = 0
    M = 1
    for i, mod in pairs:
        M *= mod
    for a, mod in pairs:
        m_i = M // mod
        result += a * inverse_mod(m_i, mod) * m_i
    return result % M

def sph_logarithm(base, target, mod):
    order = mod - 1
    prime_powers = dissolve(order)
    equations = []


    for p, exp in prime_powers.items():
        pe = p ** exp
        lookup = {power_mod(base, (order // pe) * j, mod): j for j in range(pe)}
        gamma = power_mod(target, order // pe, mod)
        for i in range(exp):
            gamma_i = power_mod(gamma, power_mod(p, i, pe), mod)
            if gamma_i in lookup:
                x_i = lookup[gamma_i]
                equations.append((x_i * power_mod(p, i, pe), pe))
                break
    return chinese_remainder_theorem(equations)


def verify_discrete_log(p, group, alpha, target):
    for x in group:
        if pow(alpha, x, p) == target:
            print(f"{alpha}^{x} ≡ {power_mod(alpha, x, p)} (mod {p})")
            print(f"Знайдено x: {x}")
            return x
    print("Значення x не знайдено.")
    return None