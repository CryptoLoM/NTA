from sph import *


if __name__ == "__main__":
    alpha = int(input("Введіть значення α: "))
    beta = int(input("Введіть значення β: "))
    p = int(input("Введіть просте число p: "))

    group = create_cyclic_group(p)
    verify = verify_discrete_log(p, group, alpha, beta)

    x_log = sph_logarithm(alpha, beta, p)
    actual_result = pow(alpha, x_log, p)
    print(f"{alpha}^{x_log} ≡ {actual_result} (mod {p})")

    start_sph = time.perf_counter()
    start_verify = time.perf_counter()
    sph_time = time.perf_counter() - start_sph
    bruteforce_time = time.perf_counter() - start_verify

    print("Час алгоритму Сільвера-Поліга-Геллмана:", sph_time)
    print("Час повного перебору:", bruteforce_time)