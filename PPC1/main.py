# PPC-1 Cálculo Numérico Aplicado. Ministrada pelo Prof. Rafael Gabler
# Mét. RK4 para resolver problema de sedimentação de uma esfera em baixo Reynolds.

import os
import math
import matplotlib.pyplot as plt

# --- Entradas fixas ----------------------
RESULTS_DIR = "results"

v0 = 0      # valores-padrão
t0 = 0


def func(x, y, st, res) -> float:
    # EDO adimensional:
    # dv/dt = (1/St) * (1 - v - (3/8) * Re * v^2)
    return (1.0 / st) * (1.0 - y - (3.0 / 8.0) * res * y ** 2)


def ks(passo, x, y, st, res) -> tuple:
    # Calcula os Ks do RK4
    k1 = func(x, y, st, res)
    k2 = func(x + 0.5 * passo, y + 0.5 * k1 * passo, st, res)
    k3 = func(x + 0.5 * passo, y + 0.5 * k2 * passo, st, res)
    k4 = func(x + passo, y + k3 * passo, st, res)
    return k1, k2, k3, k4


def rk4(passo, y, vrs: tuple) -> float:
    # Calcula o próximo y com o RK4
    k1, k2, k3, k4 = vrs
    y1 = y + (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4) * passo
    return y1


def compute_n_steps(t0: float, tf: float, h: float) -> int:
    if h <= 0:
        raise ValueError("h deve ser positivo")

    q = (tf - t0) / h
    n = round(q)

    if not math.isclose(q, n, rel_tol=0.0, abs_tol=1e-12):
        raise ValueError("h não divide bem o intervalo [t0, tf]")

    return n


def simulate(st, res, v0, t0, tf, h):
    # Roda a simulação e retorna uma lista de pares (t, v)
    n = compute_n_steps(t0, tf, h)

    t = t0
    v = v0
    results = [(t, v)]

    for i in range(1, n + 1):
        nt = t0 + h * i
        v = rk4(h, v, ks(h, t, v, st, res))
        results.append((nt, v))
        t = nt

    return results


# --- Soluções analíticas -------------------
def analytic_re0(x: float, st: float) -> tuple[float, float]:
    """
    Solução analítica para Re -> 0
        v*(t) = 1 - exp(-t / St)
    """
    return x, 1.0 - math.exp(-x / st)


def analytic_inertia(x: float, st: float, res: float) -> tuple[float, float]:
    """Solução analítica para o caso com inércia"""
    sup = -1.0 + math.sqrt(1.0 + 1.5 * res)
    inf = 0.75 * res
    v_s = sup / inf

    q = 3.0 * res / (8.0 * st)
    p = -(0.75 * res * v_s / st) - (1.0 / st)

    v = v_s + (q / p + (-1.0 / v_s - q / p) * math.exp(-p * x)) ** (-1)

    return x, v


# --- Helpers --------------------------------
def split_xy(points: list[tuple[float, float]]) -> tuple[list[float], list[float]]:
    x = [p[0] for p in points]
    y = [p[1] for p in points]
    return x, y


def max_error(numerical: list[tuple[float, float]], analytical: list[tuple[float, float]]) -> float:
    _, y_num = split_xy(numerical)
    _, y_ana = split_xy(analytical)
    return max(abs(a - b) for a, b in zip(y_num, y_ana))


def save_plot(filename: str) -> None:
    os.makedirs(RESULTS_DIR, exist_ok=True)
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, filename), dpi=300)
    plt.close()


# ITEM 1: Re -> 0: comparar numérica vs analítica para diferentes St
def item1_compare_re0():
    st_values = [0.1, 1.0, 10.0]
    res = 0.0
    h_local = 0.1
    tf_local = 60.0   # assim todos aparecem bem, inclusive St=10

    plt.figure(figsize=(9, 6))

    for st in st_values:
        numerical = simulate(st, res, v0, t0, tf_local, h_local)
        times, _ = split_xy(numerical)
        analytical = [analytic_re0(t, st) for t in times]

        err = max_error(numerical, analytical)
        print(f"[Item 1] St={st}, erro máx. = {err:.6e}")

        t_num, v_num = split_xy(numerical)
        t_ana, v_ana = split_xy(analytical)

        plt.plot(t_num, v_num, label=f"RK4, St={st}")
        plt.plot(t_ana, v_ana, "--", label=f"Analítica, St={st}")

    plt.xlabel("t*")
    plt.ylabel("v*")
    plt.title("Item 1 - Re=0: comparação para diferentes valores de St")
    plt.grid(True)
    plt.legend()
    save_plot("item1_st_combined.png")


# ITEM 2: Variar h e analisar efeito do refinamento temporal
def item2_refinement():
    st = 1.0
    res = 0.0
    tf_local = 10.0
    hs = [5.0, 2.0, 0.5, 0.1, 0.01]

    plt.figure(figsize=(9, 6))

    # curva analítica de referência, mais densa
    h_ref = 0.01
    n_ref = compute_n_steps(t0, tf_local, h_ref)
    t_exact = [t0 + i * h_ref for i in range(n_ref + 1)]
    v_exact = [analytic_re0(t, st)[1] for t in t_exact]
    plt.plot(t_exact, v_exact, "k--", linewidth=2, label="Analítica")

    for h_local in hs:
        numerical = simulate(st, res, v0, t0, tf_local, h_local)
        err = max_error(numerical, [analytic_re0(t, st) for t, _ in numerical])
        print(f"[Item 2] h={h_local}, erro máx. = {err:.6e}")

        t_num, v_num = split_xy(numerical)
        plt.plot(t_num, v_num, label=f"RK4, h={h_local}")

    plt.xlabel("t*")
    plt.ylabel("v*")
    plt.title("Item 2 - Refinamento temporal para diferentes valores de h")
    plt.grid(True)
    plt.legend()
    save_plot("item2_h_combined.png")


# ITEM 3: Resolver numericamente o caso inercial
def item3_nonzero_re():
    st = 1.0
    res = 0.5
    v0_local = 0.0
    t0_local = 0.0
    tf_local = 20.0
    h_local = 0.01

    numerical = simulate(st, res, v0_local, t0_local, tf_local, h_local)
    t_num, v_num = split_xy(numerical)

    plt.figure(figsize=(8, 5))
    plt.plot(t_num, v_num, label="Numérica (RK4)")
    plt.xlabel("t*")
    plt.ylabel("v*")
    plt.title(f"Item 3 - Solução numérica (St={st}, Re={res}, h={h_local})")
    plt.grid(True)
    plt.legend()
    save_plot("item3_re_nonzero.png")


# ITEM 4: Validar caso da inércia com solução analítica
def item4_validate_nonzero_re():
    st = 1.0
    res = 0.5
    v0_local = 0.0
    t0_local = 0.0
    tf_local = 20.0
    h_local = 0.01

    numerical = simulate(st, res, v0_local, t0_local, tf_local, h_local)
    times, _ = split_xy(numerical)
    analytical = [analytic_inertia(t, st, res) for t in times]

    err = max_error(numerical, analytical)
    print(f"[Item 4] Re={res}, erro máximo = {err:.6e}")

    t_num, v_num = split_xy(numerical)
    t_ana, v_ana = split_xy(analytical)

    plt.figure(figsize=(8, 5))
    plt.plot(t_num, v_num, label="Numérica (RK4)")
    plt.plot(t_ana, v_ana, label="Analítica")
    plt.xlabel("t*")
    plt.ylabel("v*")
    plt.title(f"Item 4 - Numérica vs Analítica (St={st}, Re={res}, h={h_local})")
    plt.grid(True)
    plt.legend()
    save_plot("item4_validation.png")


# ITEM 5: plotar para diferentes Res e discutir desvio de Re -> 0
def item5_res_sweep():
    st = 1.0
    tf_local = 20.0
    h_local = 0.1
    res_values = [0.0, 0.25, 0.5, 0.75, 1.0]

    plt.figure(figsize=(8, 5))

    for res in res_values:
        numerical = simulate(st, res, v0, t0, tf_local, h_local)
        t_num, v_num = split_xy(numerical)
        plt.plot(t_num, v_num, label=f"Re={res}")

    plt.xlabel("t*")
    plt.ylabel("v*")
    plt.title(f"Item 5 - Solução numérica para diferentes valores de Re (St={st}, h={h_local})")
    plt.grid(True)
    plt.legend()
    save_plot("item5_res_sweep.png")


def main():
    item1_compare_re0()
    item2_refinement()
    item3_nonzero_re()
    item4_validate_nonzero_re()
    item5_res_sweep()

    print("Pronto! Agora verifique a pasta 'results' :)")


if __name__ == "__main__":
    main()
