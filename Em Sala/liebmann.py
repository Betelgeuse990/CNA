"""
Cálculo Numérico Aplicado
Método das Diferenças Finitas -> Método de Liebmann
09/07/26
"""
import time

import matplotlib.pyplot as plt
import numpy as np

# Aleta
L = 1.0    # Comprimento da aleta
H = 1.0    # Altura da aleta

# Ambiente
Tb = 40    # [ºC] temperatura na base da aleta
Tamb = 25  # [ºC] temperatura ambiente

# Simulação
nodes_x = 11    # Número de Nós em X
nodes_y = 11    # Número de Nós em Y

tol = 1e-8   # Tolerância para convergência
max_it = 10_000  # Número máximo de iterações

delta_x = L / (nodes_x - 1)
delta_y = H / (nodes_y - 1)

# Grid -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def liebmann(nodes_y, nodes_x):
    # Resolve a equação de Laplace pelo método de Liebmann
    
    grid = np.full((nodes_y, nodes_x), Tamb, dtype=float)
    grid[:, 0] = Tb     # Base da aleta (x = 0)
    
    converged = False
    start_time = time.perf_counter()
    
    for it in range(1, max_it + 1):
        err = 0.0
    
        for j in range(1, nodes_y - 1):
            for i in range(1, nodes_x - 1):
                old_temp = grid[j, i]
    
                new_temp = (
                    grid[j, i + 1]      # Leste
                    + grid[j, i - 1]    # Oeste
                    + grid[j + 1, i]    # Norte
                    + grid[j - 1, i]    # Sul
                ) / 4.0
    
                grid[j, i] = new_temp
    
                local_error = abs(new_temp - old_temp)
                err = max(err, local_error)
       
        if err <= tol:    # Checa convergência
            converged = True
            break

    elapsed_time = time.perf_counter() - start_time
    
    if converged:
        print(f'Malha: {nodes_y}x{nodes_x}')
        print(f'Convergência alcançada em {it} iterações.')
        print(f'Erro final: {err:.3e}')
        print(f"Tempo computacional: {elapsed_time:.6f} s\n")
    else:
        print(f'Malha: {nodes_y}x{nodes_x}')
        print(f'O método não convergiu após {max_it} iterações.')
        print(f'Erro final: {err:.3e}')
        print(f"Tempo computacional: {elapsed_time:.6f} s\n")

    return grid, it, err, elapsed_time, converged


def post_process(grid, filename_prefix="liebmann"):
    """Salva os dados e gera os gráficos solicitados."""

    nodes_y, nodes_x = grid.shape

    x = np.linspace(0.0, L, nodes_x)
    y = np.linspace(0.0, H, nodes_y)
    X, Y = np.meshgrid(x, y)

    # Arquivo com as colunas x, y e T
    output_data = np.column_stack(
        (X.ravel(), Y.ravel(), grid.ravel())
    )

    np.savetxt(
        f"{filename_prefix}_temperature.dat",
        output_data,
        header="x y T",
        fmt="%.8e",
    )

    # Mapa bidimensional de temperatura
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
        f"{filename_prefix}_temperature_map.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)

    # Curvas de nível
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
        f"{filename_prefix}_contours.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)

    # Temperatura ao longo da linha central horizontal
    center_y_index = np.argmin(np.abs(y - H / 2.0))
    centerline_temperature = grid[center_y_index, :]

    centerline_data = np.column_stack(
        (x, centerline_temperature)
    )

    np.savetxt(
        f"{filename_prefix}_centerline.dat",
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
        title=f"Temperatura na linha central (y = {y[center_y_index]:.3f})",
        xlabel="x",
        ylabel="Temperatura (°C)",
    )

    ax.grid(True, linestyle="--", alpha=0.5)

    fig.tight_layout()
    fig.savefig(
        f"{filename_prefix}_centerline.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


def mesh_study(mesh_sizes):
    """Executa e salva os resultados do estudo de refinamento."""

    study = []
    grids = {}

    print("ESTUDO DE MALHA\n")

    for size in mesh_sizes:
        grid, iterations, error, elapsed, converged = liebmann(
            size,
            size,
        )

        x = np.linspace(0.0, L, size)
        y = np.linspace(0.0, H, size)

        center_i = np.argmin(np.abs(x - L / 2.0))
        center_j = np.argmin(np.abs(y - H / 2.0))
        center_temperature = grid[center_j, center_i]

        study.append(
            [
                size,
                size,
                center_temperature,
                iterations,
                error,
                elapsed,
                int(converged),
            ]
        )

        grids[size] = grid

    study = np.asarray(study)

    np.savetxt(
        "mesh_study.dat",
        study,
        header=(
            "nodes_x nodes_y center_temperature "
            "iterations final_error time_s converged"
        ),
        fmt=["%d", "%d", "%.8f", "%d", "%.8e", "%.8e", "%d"],
    )

    print("Resumo do estudo de malha:")
    print(
        f"{'Malha':>10} "
        f"{'T centro (°C)':>16} "
        f"{'Iterações':>12} "
        f"{'Tempo (s)':>12}"
    )

    for row in study:
        nx, ny, center_temp, iterations, _, elapsed, _ = row

        print(
            f"{int(ny)}x{int(nx):<6} "
            f"{center_temp:>16.8f} "
            f"{int(iterations):>12} "
            f"{elapsed:>12.6f}"
        )

    return grids


if __name__ == "__main__":
    mesh_sizes = [11, 21, 41, 81]

    grids = mesh_study(mesh_sizes)

    # Pós-processamento da malha mais refinada
    finest_size = mesh_sizes[-1]

    post_process(
        grids[finest_size],
        filename_prefix=f"liebmann_{finest_size}x{finest_size}",
    )
