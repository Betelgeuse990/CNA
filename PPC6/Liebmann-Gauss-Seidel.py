"""
PPC #6 — Condução bidimensional de calor em uma aleta

Equação governante:
    d²T/dx² + d²T/dy² = 0

Condições de contorno:
    x = 0: temperatura prescrita, T = Tb
    y = 0: convecção
    y = H: convecção
    x = L: convecção

Métodos:
    1. Eliminação de Gauss
    2. Liebmann / Gauss-Seidel
    3. Liebmann com relaxação
"""
import time
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# Aleta
L = 1.0       # Comprimento da aleta
H = 0.1       # Espessura da aleta
k = 200.0     # Condutividade térmica

# Ambiente
h = 25.0      # Coeficiente convectivo
Tb = 100.0    # Temperatura da base
Tamb = 25.0   # Temperatura ambiente

# Simulação
nodes_x = 11       # Número de nós em x
nodes_y = 11       # Número de nós em y

tol = 1e-8          # Tolerância para convergência
omega = 1.5         # Fator de relaxação
max_it = 100_000    # Número máximo de iterações

delta_x = L / (nodes_x - 1)
delta_y = H / (nodes_y - 1)

x = np.linspace(0.0, L, nodes_x)
y = np.linspace(0.0, H, nodes_y)


# Diretórios -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
diretorio_codigo = Path(__file__).resolve().parent

diretorio_resultados = diretorio_codigo / "resultados"

diretorio_relaxacao = (
    diretorio_resultados / "estudo_relaxacao"
)

diretorio_malha = (
    diretorio_resultados / "estudo_malha"
)

diretorio_temperatura = (
    diretorio_resultados / "campo_temperatura"
)

diretorio_linha_central = (
    diretorio_resultados / "linha_central"
)


def criar_diretorios():
    # Cria os diretórios para os resultados

    diretorios = [
        diretorio_relaxacao,
        diretorio_malha,
        diretorio_temperatura,
        diretorio_linha_central,
    ]

    for diretorio in diretorios:
        diretorio.mkdir(
            parents=True,
            exist_ok=True,
        )


# Sistema linear -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def build_system():
    # Monta a matriz A e o vetor b

    number_nodes = nodes_x * nodes_y

    A = np.zeros((number_nodes, number_nodes))
    b = np.zeros(number_nodes)

    beta = h / k

    for j in range(nodes_y):
        for i in range(nodes_x):
            p = j * nodes_x + i

            # Base da aleta
            if i == 0:
                A[p, p] = 1.0
                b[p] = Tb
                continue

            diagonal = (
                2.0 / delta_x**2
                + 2.0 / delta_y**2
            )

            # Direção x
            if i == nodes_x - 1:
                # Extremidade convectiva
                diagonal += 2.0 * beta / delta_x

                A[p, p - 1] = -2.0 / delta_x**2

                b[p] += (
                    2.0
                    * beta
                    * Tamb
                    / delta_x
                )

            else:
                A[p, p - 1] = -1.0 / delta_x**2
                A[p, p + 1] = -1.0 / delta_x**2

            # Direção y
            if j == 0:
                # Superfície inferior convectiva
                diagonal += 2.0 * beta / delta_y

                A[p, p + nodes_x] = -2.0 / delta_y**2

                b[p] += (
                    2.0
                    * beta
                    * Tamb
                    / delta_y
                )

            elif j == nodes_y - 1:
                # Superfície superior convectiva
                diagonal += 2.0 * beta / delta_y

                A[p, p - nodes_x] = -2.0 / delta_y**2

                b[p] += (
                    2.0
                    * beta
                    * Tamb
                    / delta_y
                )

            else:
                A[p, p - nodes_x] = -1.0 / delta_y**2
                A[p, p + nodes_x] = -1.0 / delta_y**2

            A[p, p] = diagonal

    return A, b


# Eliminação de Gauss -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def gaussian_elimination(A, b):
    # Resolve o sistema por eliminação de Gauss

    original_A = A.copy()
    original_b = b.copy()

    A = A.copy()
    b = b.copy()

    n = len(b)

    start_time = time.perf_counter()

    for k_index in range(n - 1):
        pivot = (
            k_index
            + np.argmax(np.abs(A[k_index:, k_index]))
        )

        if abs(A[pivot, k_index]) < 1e-14:
            raise ValueError("O sistema possui pivô nulo.")

        if pivot != k_index:
            A[[k_index, pivot]] = A[[pivot, k_index]]
            b[[k_index, pivot]] = b[[pivot, k_index]]

        for i in range(k_index + 1, n):
            factor = A[i, k_index] / A[k_index, k_index]

            A[i, k_index:] -= factor * A[k_index, k_index:]
            b[i] -= factor * b[k_index]

    temperature = np.zeros(n)

    for i in range(n - 1, -1, -1):
        known_terms = np.dot(
            A[i, i + 1:],
            temperature[i + 1:],
        )

        temperature[i] = (
            b[i] - known_terms
        ) / A[i, i]

    elapsed_time = time.perf_counter() - start_time

    error = np.max(
        np.abs(original_A @ temperature - original_b)
    )

    return {
        "temperature": temperature,
        "iterations": None,
        "error": error,
        "time": elapsed_time,
        "converged": True,
    }


# Resíduo relativo -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def relative_residual(A, b, temperature):
    # Calcula o resíduo relativo do sistema

    temperature = np.asarray(
        temperature
    ).ravel()

    numerator = np.linalg.norm(
        A @ temperature - b,
        ord=np.inf,
    )

    denominator = (
        np.linalg.norm(A, ord=np.inf)
        * np.linalg.norm(temperature, ord=np.inf)
        + np.linalg.norm(b, ord=np.inf)
    )

    return numerator / denominator


# Liebmann -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def liebmann(nodes_x, nodes_y, relaxation_factor):
    # Resolve a equação pelo método de Liebmann

    delta_x = L / (nodes_x - 1)
    delta_y = H / (nodes_y - 1)

    beta = h / k

    grid = np.full(
        (nodes_y, nodes_x),
        Tamb,
        dtype=float,
    )

    grid[:, 0] = Tb

    converged = False
    start_time = time.perf_counter()

    for it in range(1, max_it + 1):
        error = 0.0

        for j in range(nodes_y):
            for i in range(1, nodes_x):
                old_temp = grid[j, i]

                denominator = (
                    2.0 / delta_x**2
                    + 2.0 / delta_y**2
                )

                # Direção x
                if i == nodes_x - 1:
                    numerator = (
                        2.0 * grid[j, i - 1] / delta_x**2
                        + 2.0 * beta * Tamb / delta_x
                    )

                    denominator += 2.0 * beta / delta_x

                else:
                    numerator = (
                        grid[j, i - 1]
                        + grid[j, i + 1]
                    ) / delta_x**2

                # Direção y
                if j == 0:
                    numerator += (
                        2.0 * grid[j + 1, i] / delta_y**2
                        + 2.0 * beta * Tamb / delta_y
                    )

                    denominator += 2.0 * beta / delta_y

                elif j == nodes_y - 1:
                    numerator += (
                        2.0 * grid[j - 1, i] / delta_y**2
                        + 2.0 * beta * Tamb / delta_y
                    )

                    denominator += 2.0 * beta / delta_y

                else:
                    numerator += (
                        grid[j - 1, i]
                        + grid[j + 1, i]
                    ) / delta_y**2

                calculated_temp = numerator / denominator

                new_temp = (
                    relaxation_factor * calculated_temp
                    + (1.0 - relaxation_factor) * old_temp
                )

                grid[j, i] = new_temp

                local_error = abs(new_temp - old_temp)
                error = max(error, local_error)

        if error <= tol:
            converged = True
            break

    elapsed_time = time.perf_counter() - start_time

    return {
        "temperature": grid,
        "iterations": it,
        "error": error,
        "time": elapsed_time,
        "converged": converged,
    }


# Execução dos métodos -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def solve_methods(A, b):
    # Executa os três métodos

    results = {
        "gauss": gaussian_elimination(
            A,
            b,
        ),

        "liebmann": liebmann(
            nodes_x,
            nodes_y,
            relaxation_factor=1.0,
        ),

        "relaxed": liebmann(
            nodes_x,
            nodes_y,
            relaxation_factor=omega,
        ),
    }

    for result in results.values():
        result["residual"] = relative_residual(
            A,
            b,
            result["temperature"],
        )

    return results


# Resultados -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
def print_results(results):
    # Imprime os resultados e compara as soluções

    gauss = results["gauss"]
    liebmann_result = results["liebmann"]
    relaxed = results["relaxed"]

    gauss_temperature = gauss["temperature"]

    liebmann_temperature = (
        liebmann_result["temperature"].ravel()
    )

    relaxed_temperature = (
        relaxed["temperature"].ravel()
    )

    print("\nCOMPARAÇÃO DOS MÉTODOS\n")

    print(f"Malha: {nodes_y}x{nodes_x}")
    print(f"Delta x: {delta_x:.6e}")
    print(f"Delta y: {delta_y:.6e}\n")

    print(
        f"{'Método':<30}"
        f"{'Iterações':>12}"
        f"{'Resíduo relativo':>18}"
        f"{'Tempo (s)':>15}"
    )

    print("-" * 75)

    print(
        f"{'Eliminação de Gauss':<30}"
        f"{'-':>12}"
        f"{gauss['residual']:>18.3e}"
        f"{gauss['time']:>15.6f}"
    )

    print(
        f"{'Liebmann':<30}"
        f"{liebmann_result['iterations']:>12}"
        f"{liebmann_result['residual']:>18.3e}"
        f"{liebmann_result['time']:>15.6f}"
    )

    method_name = (
        f"Liebmann com omega = {omega:.2f}"
    )

    print(
        f"{method_name:<30}"
        f"{relaxed['iterations']:>12}"
        f"{relaxed['residual']:>18.3e}"
        f"{relaxed['time']:>15.6f}"
    )

    if not liebmann_result["converged"]:
        print(
            "\nLiebmann sem relaxação "
            "não convergiu."
        )

    if not relaxed["converged"]:
        print(
            "\nLiebmann com relaxação "
            "não convergiu."
        )

    difference_gauss_liebmann = np.max(
        np.abs(
            gauss_temperature
            - liebmann_temperature
        )
    )

    difference_gauss_relaxed = np.max(
        np.abs(
            gauss_temperature
            - relaxed_temperature
        )
    )

    difference_liebmann_relaxed = np.max(
        np.abs(
            liebmann_temperature
            - relaxed_temperature
        )
    )

    print("\nDIFERENÇAS ENTRE AS SOLUÇÕES\n")

    print(
        "Gauss e Liebmann: "
        f"{difference_gauss_liebmann:.3e}"
    )

    print(
        "Gauss e Liebmann relaxado: "
        f"{difference_gauss_relaxed:.3e}"
    )

    print(
        "Liebmann e Liebmann relaxado: "
        f"{difference_liebmann_relaxed:.3e}"
    )


# Solução analítica -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def analytical_solution(x):
    # Calcula a solução analítica unidimensional

    Ac = H      # Área para largura unitária
    P = 2.0     # Perímetro exposto para largura unitária

    m = np.sqrt(h * P / (k * Ac))

    numerator = (
        np.cosh(m * (L - x))
        + h / (m * k) * np.sinh(m * (L - x))
    )

    denominator = (
        np.cosh(m * L)
        + h / (m * k) * np.sinh(m * L)
    )

    theta = numerator / denominator

    return Tamb + theta * (Tb - Tamb)


# Estudo de relaxação -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
def relaxation_study(relaxation_factors):
    # Avalia diferentes fatores de relaxação

    study = []

    print("\nESTUDO DO FATOR DE RELAXAÇÃO\n")

    print(
        f"{'Omega':>10}"
        f"{'Iterações':>14}"
        f"{'Erro final':>16}"
        f"{'Tempo (s)':>14}"
    )

    print("-" * 54)

    for factor in relaxation_factors:
        result = liebmann(
            nodes_x,
            nodes_y,
            factor,
        )

        study.append(
            [
                factor,
                result["iterations"],
                result["error"],
                result["time"],
                int(result["converged"]),
            ]
        )

        print(
            f"{factor:>10.3f}"
            f"{result['iterations']:>14}"
            f"{result['error']:>16.3e}"
            f"{result['time']:>14.6f}"
        )

    study = np.asarray(study)

    np.savetxt(
        diretorio_relaxacao / "estudo_relaxacao.dat",
        study,
        header=(
            "omega iteracoes erro_final "
            "tempo_s convergiu"
        ),
        fmt=[
            "%.4f",
            "%d",
            "%.8e",
            "%.8e",
            "%d",
        ],
    )

    converged_data = study[study[:, 4] == 1]

    best_index = np.argmin(converged_data[:, 1])
    best_omega = converged_data[best_index, 0]

    print(
        "\nMelhor fator entre os valores testados: "
        f"{best_omega:.3f}"
    )

    plt.figure(figsize=(8, 5))

    plt.plot(
        converged_data[:, 0],
        converged_data[:, 1],
        marker="o",
    )

    plt.xlabel("Fator de relaxação")
    plt.ylabel("Número de iterações")
    plt.title("Efeito do fator de relaxação")
    plt.grid(True)

    plt.tight_layout()
    plt.savefig(
        diretorio_relaxacao / "estudo_relaxacao.png",
        dpi=300,
    )

    plt.close()

    return best_omega


# Estudo de malha -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
def mesh_study(mesh_sizes, relaxation_factor):
    # Avalia diferentes refinamentos de malha

    study = []

    print("\nESTUDO DE REFINAMENTO DE MALHA\n")

    print(
        f"{'Malha':>12}"
        f"{'Iterações':>14}"
        f"{'Erro médio (%)':>18}"
        f"{'Tempo (s)':>14}"
    )

    print("-" * 58)

    for size in mesh_sizes:
        result = liebmann(
            size,
            size,
            relaxation_factor,
        )

        grid = result["temperature"]

        x = np.linspace(0.0, L, size)
        y = np.linspace(0.0, H, size)

        center_j = np.argmin(
            np.abs(y - H / 2.0)
        )

        numerical_temp = grid[center_j, :]
        analytical_temp = analytical_solution(x)

        mean_error = np.mean(
            np.abs(
                (
                    numerical_temp
                    - analytical_temp
                ) / analytical_temp
            )
        ) * 100.0

        study.append(
            [
                size,
                size,
                size * size,
                result["iterations"],
                result["error"],
                mean_error,
                result["time"],
            ]
        )

        mesh_name = f"{size}x{size}"

        print(
            f"{mesh_name:>12}"
            f"{result['iterations']:>14}"
            f"{mean_error:>18.6f}"
            f"{result['time']:>14.6f}"
        )

    study = np.asarray(study)

    np.savetxt(
        diretorio_malha / "estudo_malha.dat",
        study,
        header=(
            "nos_x nos_y total_nos iteracoes "
            "erro_final erro_percentual_medio tempo_s"
        ),
        fmt=[
            "%d",
            "%d",
            "%d",
            "%d",
            "%.8e",
            "%.8e",
            "%.8e",
        ],
    )

    plt.figure(figsize=(8, 5))

    plt.plot(
        study[:, 2],
        study[:, 5],
        marker="o",
    )

    plt.xlabel("Número total de nós")
    plt.ylabel("Erro percentual médio (%)")
    plt.title("Estudo de refinamento de malha")
    plt.grid(True)

    plt.tight_layout()
    plt.savefig(
        diretorio_malha / "estudo_malha.png",
        dpi=300,
    )

    plt.close()


# Pós-processamento -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
def post_process(grid):
    # Salva os dados e gera os gráficos

    nodes_y, nodes_x = grid.shape

    x = np.linspace(0.0, L, nodes_x)
    y = np.linspace(0.0, H, nodes_y)

    X, Y = np.meshgrid(x, y)

    data = np.column_stack(
        (
            X.ravel(),
            Y.ravel(),
            grid.ravel(),
        )
    )

    np.savetxt(
        diretorio_temperatura / "campo_temperatura.dat",
        data,
        header="x y T",
        fmt="%.8e",
    )

    # Mapa de temperatura
    plt.figure(figsize=(9, 4))

    temperature_map = plt.pcolormesh(
        X,
        Y,
        grid,
        shading="auto",
        cmap="inferno",
    )

    plt.colorbar(
        temperature_map,
        label="Temperatura (°C)",
    )

    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Campo bidimensional de temperatura")

    plt.tight_layout()
    plt.savefig(
        diretorio_temperatura / "mapa_temperatura.png",
        dpi=300,
    )

    plt.close()

    # Contornos isotérmicos
    plt.figure(figsize=(9, 4))

    contours = plt.contourf(
        X,
        Y,
        grid,
        levels=20,
        cmap="inferno",
    )

    contour_lines = plt.contour(
        X,
        Y,
        grid,
        levels=10,
        colors="black",
        linewidths=0.6,
    )

    plt.clabel(
        contour_lines,
        inline=True,
        fontsize=8,
        fmt="%.1f",
    )

    plt.colorbar(
        contours,
        label="Temperatura (°C)",
    )

    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Contornos isotérmicos")

    plt.tight_layout()
    plt.savefig(
        diretorio_temperatura / "contornos_isotermicos.png",
        dpi=300,
    )

    plt.close()

    # Linha central
    center_j = np.argmin(
        np.abs(y - H / 2.0)
    )

    numerical_temp = grid[center_j, :]
    analytical_temp = analytical_solution(x)

    percent_error = np.abs(
        (
            numerical_temp
            - analytical_temp
        ) / analytical_temp
    ) * 100.0

    mean_error = np.mean(percent_error)

    centerline_data = np.column_stack(
        (
            x,
            numerical_temp,
            analytical_temp,
            percent_error,
        )
    )

    np.savetxt(
        diretorio_linha_central
    / "temperatura_linha_central.dat",
        centerline_data,
        header=(
            "x temperatura_numerica "
            "temperatura_analitica erro_percentual"
        ),
        fmt="%.8e",
    )

    plt.figure(figsize=(8, 5))

    plt.plot(
        x,
        numerical_temp,
        marker="o",
        markersize=4,
        label="Solução numérica",
    )

    plt.plot(
        x,
        analytical_temp,
        linestyle="--",
        label="Solução analítica",
    )

    plt.xlabel("x")
    plt.ylabel("Temperatura (°C)")
    plt.title("Temperatura ao longo da linha central")
    plt.grid(True)
    plt.legend()

    plt.tight_layout()
    plt.savefig(
        diretorio_linha_central
    / "temperatura_linha_central.png",
        dpi=300,
    )

    plt.close()

    print(
        f"\nErro percentual médio na linha central: "
        f"{mean_error:.6f}%"
    )


"""
As discussões do final do PPC em si estarão no README do github.
"""
if __name__ == "__main__":
    criar_diretorios()

    A, b = build_system()

    results = solve_methods(
        A,
        b,
    )

    print_results(results)

    relaxation_factors = [
        0.80,
        1.00,
        1.20,
        1.40,
        1.60,
        1.80,
        1.90,
        1.94,
        1.95,
        1.955,
        1.96,
        1.965,
        1.97,
        1.98,
    ]

    best_omega = relaxation_study(
        relaxation_factors
    )

    mesh_sizes = [
        11,
        21,
        31,
    ]

    mesh_study(
        mesh_sizes,
        best_omega,
    )

    final_size = mesh_sizes[-1]

    final_result = liebmann(
        final_size,
        final_size,
        best_omega,
    )

    post_process(
        final_result["temperature"]
    )

    print(
        "\nResultados salvos em:"
        f"\n{diretorio_resultados}"
    )
