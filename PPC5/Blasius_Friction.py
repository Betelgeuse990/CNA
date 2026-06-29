"""
PPC #5 - Resolução Numérica da Equação de Blasius

Métodos utilizados:
1. Método do Tiro para transformar o PVC em PVI;
2. Método da Secante para atualizar o parâmetro de tiro;
3. Runge-Kutta de quarta ordem para integrar o sistema de EDOs.
"""

import math
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt


# ============================================================================
# SISTEMA DE EQUAÇÕES DE BLASIUS
# ============================================================================

def f1(y: np.ndarray) -> float:
    """
    F1(y1, y2, y3) = dy1/deta.
    """
    return y[1]


def f2(y: np.ndarray) -> float:
    """
    F2(y1, y2, y3) = dy2/deta.
    """
    return y[2]


def f3(y: np.ndarray) -> float:
    """
    F3(y1, y2, y3) = dy3/deta.
    """
    return -0.5 * y[0] * y[2]


def sistema_blasius(y: np.ndarray) -> np.ndarray:
    """
    Calcula simultaneamente as derivadas do sistema de Blasius.

    Parâmetros
    ----------
    y : np.ndarray
        Vetor [y1, y2, y3], equivalente a [f, f', f''].

    Retorna
    -------
    np.ndarray
        Vetor [dy1/deta, dy2/deta, dy3/deta].
    """
    return np.array([
        f1(y),
        f2(y),
        f3(y),
    ], dtype=float)


# ============================================================================
# RUNGE-KUTTA DE QUARTA ORDEM
# ============================================================================

def coefs_rk4(
    step: float,
    y: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Calcula os quatro vetores de coeficientes do método RK4.

    Convenção adotada:
        K1 = h * F(y_i)
        K2 = h * F(y_i + K1 / 2)
        K3 = h * F(y_i + K2 / 2)
        K4 = h * F(y_i + K3)

    Os vetores K já incorporam o passo de integração h.
    """
    y = np.asarray(y, dtype=float)

    K1 = step * sistema_blasius(y)

    K2 = step * sistema_blasius(
        y + K1 / 2.0
    )

    K3 = step * sistema_blasius(
        y + K2 / 2.0
    )

    K4 = step * sistema_blasius(
        y + K3
    )

    return K1, K2, K3, K4


def step_rk4(
    step: float,
    y: np.ndarray,
) -> np.ndarray:
    """
    Executa um passo completo do método RK4 para o sistema de Blasius.

    Retorna o vetor atualizado [y1, y2, y3].
    """
    K1, K2, K3, K4 = coefs_rk4(step, y)

    y_novo = y + (
        K1
        + 2.0 * K2
        + 2.0 * K3
        + K4
    ) / 6.0

    return y_novo


# ============================================================================
# INTEGRAÇÃO DO PVI PARA UM DADO CHUTE s = f''(0)
# ============================================================================

def integrate_rk4(
    step: float,
    eta_max: float,
    s: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Integra o sistema de Blasius de eta = 0 até eta = eta_max.

    Parâmetros
    ----------
    step : float
        Passo de integração h = delta_eta.
    eta_max : float
        Limite superior da variável de similaridade.
    s : float
        Parâmetro de tiro, correspondente a y3(0) = f''(0).

    Retorna
    -------
    eta_values : np.ndarray
        Vetor com os valores de eta.
    y_values : np.ndarray
        Matriz em que cada linha contém [y1, y2, y3],
        ou seja, [f, f', f''] para cada eta.
    """
    if step <= 0:
        raise ValueError("O passo de integração deve ser positivo.")

    if eta_max <= 0:
        raise ValueError("eta_max deve ser positivo.")

    eta = 0.0

    # Condições iniciais do PVI:
    # y1(0) = f(0) = 0
    # y2(0) = f'(0) = 0
    # y3(0) = f''(0) = s
    y = np.array([0.0, 0.0, s], dtype=float)

    eta_values = [eta]
    y_values = [y.copy()]

    while eta < eta_max:
        # Evita ultrapassar eta_max no último passo.
        h_atual = min(step, eta_max - eta)

        y = step_rk4(h_atual, y)
        eta += h_atual

        eta_values.append(eta)
        y_values.append(y.copy())

    return np.array(eta_values), np.array(y_values)


# ============================================================================
# AVALIAÇÃO DE UM CHUTE DO MÉTODO DO TIRO
# ============================================================================

def avaliar_chute(
    s: float,
    delta_eta: float,
    eta_max: float,
) -> tuple[float, np.ndarray, np.ndarray]:
    """
    Integra o sistema para um chute s e calcula o erro do Método do Tiro.

    A função erro é:
        E(s) = f'(eta_max) - 1

    Retorna
    -------
    erro : float
        Erro associado ao chute avaliado.
    eta_values : np.ndarray
        Vetor de eta.
    y_values : np.ndarray
        Perfil completo [f, f', f''].
    """
    eta_values, y_values = integrate_rk4(
        delta_eta,
        eta_max,
        s,
    )

    erro = y_values[-1, 1] - 1.0

    return erro, eta_values, y_values


# ============================================================================
# MÉTODO DO TIRO COM ATUALIZAÇÃO PELA SECANTE
# ============================================================================

def metodo_tiro_secante(
    s0: float,
    s1: float,
    delta_eta: float,
    eta_max: float,
    tol: float,
    max_tiro: int,
) -> tuple[float, float, int, np.ndarray, np.ndarray]:
    """
    Resolve a equação de Blasius pelo Método do Tiro.

    O parâmetro de tiro s = f''(0) é atualizado pelo método da secante.

    Retorna
    -------
    s_convergido : float
        Valor numérico encontrado para f''(0).
    erro_final : float
        E(s) = f'(eta_max) - 1 para o último chute.
    iteracoes : int
        Número de atualizações da secante realizadas.
    eta_values : np.ndarray
        Valores de eta associados ao perfil convergido.
    y_values : np.ndarray
        Perfil convergido [f, f', f''].
    """
    erro_s0, _, _ = avaliar_chute(
        s0,
        delta_eta,
        eta_max,
    )

    erro_s1, eta_values, y_values = avaliar_chute(
        s1,
        delta_eta,
        eta_max,
    )

    # Caso o segundo chute já satisfaça a tolerância.
    if abs(erro_s1) < tol:
        return s1, erro_s1, 0, eta_values, y_values

    for it in range(1, max_tiro + 1):

        denominador = erro_s1 - erro_s0

        if np.isclose(denominador, 0.0):
            raise ZeroDivisionError(
                "A secante não pode avançar: os erros dos dois chutes "
                "são muito próximos."
            )

        # Fórmula da secante.
        s2 = s1 - erro_s1 * (s1 - s0) / denominador

        # Novo chute: reinicia a integração em eta = 0.
        erro_s2, eta_values, y_values = avaliar_chute(
            s2,
            delta_eta,
            eta_max,
        )

        if abs(erro_s2) < tol:
            return s2, erro_s2, it, eta_values, y_values

        # Atualiza os dois pares mais recentes:
        # (s0, erro_s0) e (s1, erro_s1).
        s0 = s1
        erro_s0 = erro_s1

        s1 = s2
        erro_s1 = erro_s2

    raise RuntimeError(
        "O Método do Tiro não convergiu dentro do número máximo "
        "de iterações."
    )


# ============================================================================
# PÓS-PROCESSAMENTO: eta_99 E ATRITO NA PAREDE
# ============================================================================

def calcular_eta99(
    eta_values: np.ndarray,
    y_values: np.ndarray,
) -> float:
    """
    Determina eta_99 por interpolação linear.

    A definição utilizada é:
        f'(eta_99) = 0.99
    """
    f_linha = y_values[:, 1]
    alvo = 0.99

    indices = np.where(
        (f_linha[:-1] <= alvo)
        & (f_linha[1:] >= alvo)
    )[0]

    if len(indices) == 0:
        raise ValueError(
            "O perfil não atingiu f'(eta) = 0.99. "
            "Aumente eta_max."
        )

    i = indices[0]

    eta_a = eta_values[i]
    eta_b = eta_values[i + 1]

    f_a = f_linha[i]
    f_b = f_linha[i + 1]

    # Caso raro em que um ponto calculado seja exatamente 0.99.
    if np.isclose(f_a, alvo):
        return eta_a

    eta99 = eta_a + (
        (alvo - f_a)
        * (eta_b - eta_a)
        / (f_b - f_a)
    )

    return eta99


def calcular_atrito(
    s_convergido: float,
    u_inf: float,
    x: float,
    nu: float,
) -> tuple[float, float]:
    """
    Calcula o número de Reynolds local e o coeficiente de atrito local.

    Re_x = U_inf * x / nu

    C_f = 2 * f''(0) / sqrt(Re_x)
    """
    if u_inf <= 0 or x <= 0 or nu <= 0:
        raise ValueError(
            "U_inf, x e nu devem ser positivos."
        )

    re_x = u_inf * x / nu

    cf = 2.0 * s_convergido / math.sqrt(re_x)

    return re_x, cf


# ============================================================================
# SALVAMENTO DOS RESULTADOS E GERAÇÃO DE GRÁFICOS
# ============================================================================

def salvar_perfil_csv(
    pasta_resultados: Path,
    eta_values: np.ndarray,
    y_values: np.ndarray,
) -> Path:
    """
    Salva o perfil numérico da solução em arquivo CSV.

    O arquivo contém as colunas:
        eta, f, f_linha, f_duas_linhas
    """
    arquivo_csv = pasta_resultados / "perfil_blasius.csv"

    dados = np.column_stack((
        eta_values,
        y_values[:, 0],
        y_values[:, 1],
        y_values[:, 2],
    ))

    np.savetxt(
        arquivo_csv,
        dados,
        delimiter=",",
        header="eta,f,f_linha,f_duas_linhas",
        comments="",
        fmt="%.10e",
    )

    return arquivo_csv


def gerar_graficos(
    pasta_resultados: Path,
    eta_values: np.ndarray,
    y_values: np.ndarray,
    eta99: float,
) -> list[Path]:
    """
    Gera e salva os três gráficos exigidos no PPC:

    1. f(eta);
    2. f'(eta);
    3. f''(eta).

    Retorna uma lista com os caminhos dos arquivos gerados.
    """
    arquivos_gerados = []

    # ------------------------------------------------------------------------
    # Gráfico de f(eta)
    # ------------------------------------------------------------------------
    arquivo_f = pasta_resultados / "perfil_f.png"

    plt.figure(figsize=(8, 5))
    plt.plot(eta_values, y_values[:, 0])

    plt.xlabel(r"$\eta$")
    plt.ylabel(r"$f(\eta)$")
    plt.title(r"Perfil de similaridade $f(\eta)$")

    plt.grid(True)
    plt.tight_layout()

    plt.savefig(arquivo_f, dpi=300)
    plt.close()

    arquivos_gerados.append(arquivo_f)

    # ------------------------------------------------------------------------
    # Gráfico de f'(eta)
    # ------------------------------------------------------------------------
    arquivo_f_linha = pasta_resultados / "perfil_f_linha.png"

    plt.figure(figsize=(8, 5))
    plt.plot(
        eta_values,
        y_values[:, 1],
        label=r"$f'(\eta)=u/U_\infty$",
    )

    plt.axhline(
        1.0,
        linestyle="--",
        label=r"Limite: $f'(\eta)=1$",
    )

    plt.axhline(
        0.99,
        linestyle="--",
        label=r"Critério: $f'(\eta)=0.99$",
    )

    plt.axvline(
        eta99,
        linestyle="--",
        label=rf"$\eta_{{99}}={eta99:.4f}$",
    )

    plt.xlabel(r"$\eta$")
    plt.ylabel(r"$f'(\eta)=u/U_\infty$")
    plt.title("Perfil de velocidade adimensional")

    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.savefig(arquivo_f_linha, dpi=300)
    plt.close()

    arquivos_gerados.append(arquivo_f_linha)

    # ------------------------------------------------------------------------
    # Gráfico de f''(eta)
    # ------------------------------------------------------------------------
    arquivo_f_duas_linhas = pasta_resultados / "perfil_f_duas_linhas.png"

    plt.figure(figsize=(8, 5))
    plt.plot(eta_values, y_values[:, 2])

    plt.xlabel(r"$\eta$")
    plt.ylabel(r"$f''(\eta)$")
    plt.title(r"Perfil de similaridade $f''(\eta)$")

    plt.grid(True)
    plt.tight_layout()

    plt.savefig(arquivo_f_duas_linhas, dpi=300)
    plt.close()

    arquivos_gerados.append(arquivo_f_duas_linhas)

    return arquivos_gerados


def salvar_resumo(
    pasta_resultados: Path,
    s_convergido: float,
    erro_final: float,
    iteracoes: int,
    f_linha_final: float,
    eta99: float,
    re_x: float,
    cf: float,
    delta_eta: float,
    eta_max: float,
    tol: float,
) -> Path:
    """
    Salva um resumo textual dos principais resultados numéricos.
    """
    arquivo_resumo = pasta_resultados / "resumo_resultados.txt"

    valor_referencia_s = 0.332057
    valor_referencia_c_delta = 4.92

    erro_s_referencia = abs(
        s_convergido - valor_referencia_s
    )

    erro_c_delta = abs(
        eta99 - valor_referencia_c_delta
    )

    with open(arquivo_resumo, "w", encoding="utf-8") as arquivo:
        arquivo.write("PPC #5 - Solução Numérica da Equação de Blasius\n")
        arquivo.write("=" * 58 + "\n\n")

        arquivo.write("Parâmetros numéricos\n")
        arquivo.write("-" * 58 + "\n")
        arquivo.write(f"Passo de integração delta_eta: {delta_eta:.10e}\n")
        arquivo.write(f"eta_max: {eta_max:.10f}\n")
        arquivo.write(f"Tolerância: {tol:.10e}\n\n")

        arquivo.write("Método do Tiro\n")
        arquivo.write("-" * 58 + "\n")
        arquivo.write(f"f''(0) convergido: {s_convergido:.10f}\n")
        arquivo.write(
            f"Valor de referência de f''(0): "
            f"{valor_referencia_s:.10f}\n"
        )
        arquivo.write(
            f"Erro absoluto em f''(0): "
            f"{erro_s_referencia:.10e}\n"
        )
        arquivo.write(f"Iterações da secante: {iteracoes}\n")
        arquivo.write(f"Erro final do tiro: {erro_final:.10e}\n")
        arquivo.write(
            f"f'(eta_max): {f_linha_final:.10f}\n\n"
        )

        arquivo.write("Camada limite\n")
        arquivo.write("-" * 58 + "\n")
        arquivo.write(f"eta_99: {eta99:.10f}\n")
        arquivo.write(f"C_delta = eta_99: {eta99:.10f}\n")
        arquivo.write(
            f"Referência para C_delta: "
            f"{valor_referencia_c_delta:.10f}\n"
        )
        arquivo.write(
            f"Erro absoluto em C_delta: "
            f"{erro_c_delta:.10e}\n\n"
        )

        arquivo.write("Atrito na parede\n")
        arquivo.write("-" * 58 + "\n")
        arquivo.write(f"Re_x: {re_x:.10e}\n")
        arquivo.write(f"C_f: {cf:.10e}\n")

    return arquivo_resumo


# ============================================================================
# EXECUÇÃO PRINCIPAL
# ============================================================================

def main() -> None:
    # Parâmetros do Método do Tiro
    s0 = 0.30
    s1 = 0.40

    # Parâmetros da integração numérica
    delta_eta = 1e-3
    eta_max = 10.0
    tol = 1e-8
    max_tiro = 50

    # Dados físicos para o cálculo de Re_x e C_f.
    # Exemplo: ar em condição ambiente.
    u_inf = 5.0       # m/s
    x = 1.0           # m
    nu = 1.5e-5       # m²/s

    (
        s_convergido,
        erro_final,
        iteracoes,
        eta_values,
        y_values,
    ) = metodo_tiro_secante(
        s0,
        s1,
        delta_eta,
        eta_max,
        tol,
        max_tiro,
    )

    eta99 = calcular_eta99(
        eta_values,
        y_values,
    )

    re_x, cf = calcular_atrito(
        s_convergido,
        u_inf,
        x,
        nu,
    )

    # ------------------------------------------------------------------------
    # Criação da pasta de resultados
    # ------------------------------------------------------------------------
    pasta_resultados = Path("resultados")
    pasta_resultados.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------------
    # Salvamento do perfil numérico
    # ------------------------------------------------------------------------
    arquivo_csv = salvar_perfil_csv(
        pasta_resultados,
        eta_values,
        y_values,
    )

    # ------------------------------------------------------------------------
    # Geração dos gráficos
    # ------------------------------------------------------------------------
    arquivos_graficos = gerar_graficos(
        pasta_resultados,
        eta_values,
        y_values,
        eta99,
    )

    # ------------------------------------------------------------------------
    # Salvamento do resumo final
    # ------------------------------------------------------------------------
    arquivo_resumo = salvar_resumo(
        pasta_resultados,
        s_convergido,
        erro_final,
        iteracoes,
        y_values[-1, 1],
        eta99,
        re_x,
        cf,
        delta_eta,
        eta_max,
        tol,
    )

    print("\nResultados finais")
    print("-" * 45)

    print(f"f''(0) convergido: {s_convergido:.10f}")
    print("Valor de referência: 0.3320570000")

    print(
        "Erro absoluto em relação à referência: "
        f"{abs(s_convergido - 0.332057):.10e}"
    )

    print(f"Iterações da secante: {iteracoes}")
    print(f"Erro final do tiro: {erro_final:.10e}")

    print(
        f"f'(eta_max = {eta_max:.2f}): "
        f"{y_values[-1, 1]:.10f}"
    )

    print(f"eta_99: {eta99:.10f}")
    print(f"C_delta = eta_99: {eta99:.10f}")
    print("Referência para C_delta: 4.9200000000")

    print(f"Re_x: {re_x:.6e}")
    print(f"C_f: {cf:.6e}")

    print("\nArquivos gerados")
    print("-" * 45)
    print(f"Perfil numérico: {arquivo_csv}")
    print(f"Resumo numérico: {arquivo_resumo}")

    for arquivo_grafico in arquivos_graficos:
        print(f"Gráfico: {arquivo_grafico}")

if __name__ == "__main__":
    main()
