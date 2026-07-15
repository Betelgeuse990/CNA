"""
Cálculo Numérico Aplicado
Método das Diferenças Finitas -> Método de Liebmann
09/07/26
"""

import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def liebmann(L, H, nodes_x, nodes_y, Tb, Tamb, tol, max_it):
    """Resolve a equação de Laplace 2D pelo método de Liebmann."""

    if L <= 0 or H <= 0:
        raise ValueError("L e H devem ser positivos.")

    if nodes_x < 3 or nodes_y < 3:
        raise ValueError(
            "A malha deve possuir pelo menos 3 nós em cada direção."
        )

    if tol <= 0 or max_it < 1:
        raise ValueError(
            "A tolerância e o número máximo de iterações devem ser positivos."
        )

    delta_x = L / (nodes_x - 1)
    delta_y = H / (nodes_y - 1)

    # Estimativa inicial:
    # todos os nós começam com a temperatura ambiente.
    grid = np.full(
        (nodes_y, nodes_x),
        Tamb,
        dtype=float,
    )

    # Condições de contorno de Dirichlet:
    # fronteira esquerda na temperatura da base;
    # demais fronteiras permanecem em Tamb.
    grid[:, 0] = Tb

    # Coeficientes da equação discretizada.
    inv_dx2 = 1.0 / delta_x**2
    inv_dy2 = 1.0 / delta_y**2
    denominator = 2.0 * (inv_dx2 + inv_dy2)

    converged = False
    err = np.inf

    start_time = time.perf_counter()

    # Método de Liebmann (Gauss-Seidel).
    for iteration in range(1, max_it + 1):
        err = 0.0

        for j in range(1, nodes_y - 1):
            for i in range(1, nodes_x - 1):
                old_temp = grid[j, i]

                # Forma geral da equação de Laplace discretizada.
                # Também funciona quando delta_x != delta_y.
                new_temp = (
                    (
                        grid[j, i + 1]
                        + grid[j, i - 1]
                    ) * inv_dx2
                    + (
                        grid[j + 1, i]
                        + grid[j - 1, i]
                    ) * inv_dy2
                ) / denominator

                # Atualização imediata característica de Gauss-Seidel.
                grid[j, i] = new_temp

                local_error = abs(new_temp - old_temp)
                err = max(err, local_error)

        if err <= tol:
            converged = True
            break

    elapsed_time = time.perf_counter() - start_time

    return (
        grid,
        iteration,
        err,
        elapsed_time,
        converged,
    )


def post_process(
    grid,
    L,
    H,
    output_dir,
    filename_prefix="liebmann",
):
    """Salva os resultados e gera os gráficos solicitados."""

    output_dir = Path(output_dir)
    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    nodes_y, nodes_x = grid.shape

    x = np.linspace(0.0, L, nodes_x)
    y = np.linspace(0.0, H, nodes_y)

    X, Y = np.meshgrid(x, y)

    # =====================================================
    # Arquivo com as colunas x, y e T
    # =====================================================

    output_data = np.column_stack(
        (
            X.ravel(),
            Y.ravel(),
            grid.ravel(),
        )
    )

    np.savetxt(
        output_dir / f"{filename_prefix}_temperature.dat",
        output_data,
        header="x y T",
        fmt="%.8e",
    )

    # =====================================================
    # Mapa bidimensional de temperatura
    # =====================================================

    fig, ax = plt.subplots(figsize=(8, 5))

    temperature_map = ax.pcolormesh(
        X,
        Y,
        grid,
        shading="auto",
        cmap="inferno",
    )

    fig.colorbar(
        temperature_map,
        ax=ax,
        label="Temperatura (°C)",
    )

    ax.set(
        title="Campo bidimensional de temperatura",
        xlabel="x",
        ylabel="y",
        aspect="equal",
    )

    fig.tight_layout()

    fig.savefig(
        output_dir / f"{filename_prefix}_temperature_map.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)

    # =====================================================
    # Curvas de nível
    # =====================================================

    fig, ax = plt.subplots(figsize=(8, 5))

    filled_contours = ax.contourf(
        X,
        Y,
        grid,
        levels=20,
        cmap="inferno",
    )

    contour_lines = ax.contour(
        X,
        Y,
        grid,
        levels=10,
        colors="black",
        linewidths=0.6,
    )

    ax.clabel(
        contour_lines,
        inline=True,
        fontsize=8,
        fmt="%.1f °C",
    )

    fig.colorbar(
        filled_contours,
        ax=ax,
        label="Temperatura (°C)",
    )

    ax.set(
        title="Contornos isotérmicos",
        xlabel="x",
        ylabel="y",
        aspect="equal",
    )

    fig.tight_layout()

    fig.savefig(
        output_dir / f"{filename_prefix}_contours.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)

    # =====================================================
    # Temperatura ao longo da linha central horizontal
    # =====================================================

    center_y_index = np.argmin(
        np.abs(y - H / 2.0)
    )

    centerline_temperature = grid[
        center_y_index,
        :,
    ]

    centerline_data = np.column_stack(
        (
            x,
            centerline_temperature,
        )
    )

    np.savetxt(
        output_dir / f"{filename_prefix}_centerline.dat",
        centerline_data,
        header="x T",
        fmt="%.8e",
    )

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.plot(
        x,
        centerline_temperature,
        color="tab:red",
        linewidth=2,
        marker="o",
        markersize=3,
    )

    ax.set(
        title=(
            "Temperatura na linha central "
            f"(y = {y[center_y_index]:.3f})"
        ),
        xlabel="x",
        ylabel="Temperatura (°C)",
    )

    ax.grid(
        True,
        linestyle="--",
        alpha=0.5,
    )

    fig.tight_layout()

    fig.savefig(
        output_dir / f"{filename_prefix}_centerline.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)


def mesh_study(
    mesh_sizes,
    L,
    H,
    Tb,
    Tamb,
    tol,
    max_it,
    output_dir,
):
    """Executa o estudo de refinamento de malha."""

    output_dir = Path(output_dir)
    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    study = []
    grids = {}

    print("\nESTUDO DE MALHA")

    print(
        f"{'Malha':>11} "
        f"{'T centro (°C)':>16} "
        f"{'Iterações':>12} "
        f"{'Erro final':>14} "
        f"{'Tempo (s)':>12} "
        f"{'Status':>14}"
    )

    for nodes_x, nodes_y in mesh_sizes:
        (
            grid,
            iterations,
            error,
            elapsed,
            converged,
        ) = liebmann(
            L,
            H,
            nodes_x,
            nodes_y,
            Tb,
            Tamb,
            tol,
            max_it,
        )

        x = np.linspace(
            0.0,
            L,
            nodes_x,
        )

        y = np.linspace(
            0.0,
            H,
            nodes_y,
        )

        center_i = np.argmin(
            np.abs(x - L / 2.0)
        )

        center_j = np.argmin(
            np.abs(y - H / 2.0)
        )

        center_temperature = grid[
            center_j,
            center_i,
        ]

        study.append(
            [
                nodes_x,
                nodes_y,
                center_temperature,
                iterations,
                error,
                elapsed,
                int(converged),
            ]
        )

        grids[(nodes_x, nodes_y)] = grid

        if converged:
            status = "convergiu"
        else:
            status = "não convergiu"

        mesh_name = f"{nodes_x}x{nodes_y}"

        print(
            f"{mesh_name:>11} "
            f"{center_temperature:>16.8f} "
            f"{iterations:>12d} "
            f"{error:>14.3e} "
            f"{elapsed:>12.6f} "
            f"{status:>14}"
        )

    study = np.asarray(study)

    np.savetxt(
        output_dir / "mesh_study.dat",
        study,
        header=(
            "nodes_x nodes_y center_temperature "
            "iterations final_error time_s converged"
        ),
        fmt=[
            "%d",
            "%d",
            "%.8f",
            "%d",
            "%.8e",
            "%.8e",
            "%d",
        ],
    )

    return grids, study


def read_float(
    message,
    default,
    positive=False,
):
    """
    Lê um número real.

    Pressionar Enter mantém o valor padrão.
    """

    while True:
        text = input(
            f"{message} [{default}]: "
        ).strip()

        # Permite utilizar vírgula como separador decimal.
        text = text.replace(",", ".")

        try:
            if text == "":
                value = float(default)
            else:
                value = float(text)

            if positive and value <= 0:
                raise ValueError

            return value

        except ValueError:
            if positive:
                print(
                    "Digite um número válido e positivo."
                )
            else:
                print(
                    "Digite um número válido."
                )


def read_int(
    message,
    default,
    minimum=1,
):
    """
    Lê um número inteiro.

    Pressionar Enter mantém o valor padrão.
    """

    while True:
        text = input(
            f"{message} [{default}]: "
        ).strip()

        try:
            if text == "":
                value = int(default)
            else:
                value = int(text)

            if value < minimum:
                raise ValueError

            return value

        except ValueError:
            print(
                "Digite um número inteiro "
                f"maior ou igual a {minimum}."
            )


def main():
    print(
        "CONDUÇÃO BIDIMENSIONAL EM UMA ALETA"
    )

    print(
        "Pressione Enter para aceitar os valores "
        "indicados entre colchetes.\n"
    )

    # =====================================================
    # Entrada dos dados
    # =====================================================

    L = read_float(
        "Comprimento L",
        1.0,
        positive=True,
    )

    H = read_float(
        "Altura H",
        1.0,
        positive=True,
    )

    nodes_x = read_int(
        "Número inicial de nós em x",
        11,
        minimum=3,
    )

    nodes_y = read_int(
        "Número inicial de nós em y",
        11,
        minimum=3,
    )

    Tb = read_float(
        "Temperatura da base Tb (°C)",
        40.0,
    )

    Tamb = read_float(
        "Temperatura ambiente T_inf (°C)",
        25.0,
    )

    tol = read_float(
        "Tolerância",
        1e-8,
        positive=True,
    )

    max_it = read_int(
        "Número máximo de iterações",
        50_000,
        minimum=1,
    )

    # =====================================================
    # Estudo de refinamento
    # =====================================================

    # Para uma malha inicial 11x11, produz:
    # 11x11, 21x21, 41x41 e 81x81.
    refinement_factors = (
        1,
        2,
        4,
        8,
    )

    mesh_sizes = [
        (
            1 + (nodes_x - 1) * factor,
            1 + (nodes_y - 1) * factor,
        )
        for factor in refinement_factors
    ]

    output_dir = Path(
        "resultados_liebmann"
    )

    grids, study = mesh_study(
        mesh_sizes,
        L,
        H,
        Tb,
        Tamb,
        tol,
        max_it,
        output_dir,
    )

    # =====================================================
    # Pós-processamento da malha mais refinada
    # =====================================================

    finest_nodes_x, finest_nodes_y = mesh_sizes[-1]

    finest_grid = grids[
        (
            finest_nodes_x,
            finest_nodes_y,
        )
    ]

    prefix = (
        f"liebmann_"
        f"{finest_nodes_x}x{finest_nodes_y}"
    )

    post_process(
        finest_grid,
        L,
        H,
        output_dir,
        prefix,
    )

    print(
        "\nResultados salvos em:"
    )

    print(
        output_dir.resolve()
    )


if __name__ == "__main__":
    main()
