"""
PPC #5 - Resolução Numérica da Equação de Blasius
Autor: Estevão André
Docente: Prof. Dr. Rafael Gabler

Métodos usados:
1. Tiro para transformação do problema de PVC para PVI;
2. Runge-Kutta de 4ª ordem para integração do sistema de EDOs
"""
import numpy as np

# Header
header = '''PPC #5 - Blasius Equation Solver
by Estevão A.
'''
print(header + '=-' * 30)

# Entradas -> obs: não quis deixar tudo com input() para agilizar...
s0 = 0.30
s1 = 0.40

delta_n = 1e-3
n_max = 10.0

tol = 1e-8
max_tiro = 50

# Condições iniciais
y1_0 = y2_0 = 0.0

# Funções F1, F2 e F3 do sistema de EDOs
# Recebem um vetor com os valores de y1, y2 e y3
def f1(y) -> float:
    y1 = y[0]
    y2 = y[1]
    y3 = y[2]
    return y2

def f2(y) -> float:
    y1 = y[0]
    y2 = y[1]
    y3 = y[2]
    return y3

def f3(y) -> float:
    y1 = y[0]
    y2 = y[1]
    y3 = y[2]
    return - 0.5 * y1 * y3


def sistema_blasius(y: np.ndarray) -> np.ndarray:
    """
    Calcula simultaneamente as derivadas do sistema de Blasius.

    Parâmetros
    ----------
    y : np.ndarray
        Vetor [y1, y2, y3].

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


def coefs_rk4(step: float, y: np.ndarray) -> tuple[np.ndarray, ...]:
    """
    Calcula os quatro vetores de coeficientes do método RK4.

    Convenção adotada:
        K1 = h * F(y_i)
        K2 = h * F(y_i + K1/2)
        K3 = h * F(y_i + K2/2)
        K4 = h * F(y_i + K3)

    Portanto, os coeficientes já incorporam o passo de integração.
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


def step_rk4(step: float, y: np.ndarray) -> np.ndarray:
    """
    Executa um passo completo do método de Runge-Kutta de quarta ordem
    para o sistema de três EDOs da equação de Blasius.
    """
    K1, K2, K3, K4 = coefs_rk4(step, y)

    y_novo = y + (
        K1
        + 2.0 * K2
        + 2.0 * K3
        + K4
    ) / 6.0

    return y_novo


def integrate_rk4(step, y: np.ndarray) -> np.ndarray:
    h = step
    n = 0.0     # Ponto de partida da integração

    while n < n_max:
        y = step_rk4(h, y)
        n += h

    return y


# Usando chute s0:
vec_y = np.array([y1_0, y2_0, s0])
y_s0 = integrate_rk4(delta_n, vec_y)
err_s0 = y_s0[1] - 1

# Usando chute s1:
vec_y = np.array([y1_0, y2_0, s1])
y_s1 = integrate_rk4(delta_n, vec_y)
err_s1 = y_s1[1] - 1

# M. da Secante
it = 0  # iterações da secante

while abs(err_s1) > tol and it < max_tiro:

    if err_s1 == err_s0:
        raise ZeroDivisionError(
            'A secante não pode avançar: os dois chutes são iguais.'
        )

    next_s = s1 - err_s1 * (s1 - s0) / (err_s1 - err_s0)

    s0 = s1
    err_s0 = err_s1

    s1 = next_s

    vec_y_s1 = np.array([y1_0, y2_0, s1], dtype=float)
    y_s1 = integrate_rk4(delta_n, vec_y_s1)

    err_s1 = y_s1[1] - 1.0
    it += 1

print(y_s1)
print(f'Valor de S que convergiu: {s1}')
print(f'Erro: {err_s1}')
