"""
P1 - CNA
Estudante: Estevão A.
"""

import numpy as np


err = 1e-10     # Erro máximo tolerado para algo ser considerado como 0.0

# Dados experimentais fornecidos
Q = np.array(
    [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
    )

Qt = np.array(
    [18.2, 26.8, 34.7, 41.5, 47.0, 51.2, 54.1]
    )

Pb = np.array(
    [0.9, 1.5, 2.6, 4.1, 6.2, 8.9, 12.1]
)


"""
Questão 1 - Resolvendo sistemas lineares
Implementação: matrizes LU
"""

def solver_LU(mat_A, vec_B) -> np.array:
    """
    Recebe uma matriz dos coeficientes mat_A e um vetor resposta vec_B e
    resolve: mat_A . X = vec_B e retorna o vetor-solução "x".
    """
    size = len(vec_B)

    D = np.zeros(size)
    X = np.zeros(size)

    n = 3    # gambiarra de adaptação
    _n = 3
    A = mat_A[:]
    B = vec_B[:]

    # --- Eliminação progressiva =====================
    def eliminacao(matrix):
        _n = len(matrix)
        U = np.array(matrix, dtype=float)  # A tal matriz "U" (cópia em float)
        L = np.array([[1.0 if i == j else 0.0 for j in range(_n)] for i in range(_n)])    # A tal matriz "L", identidade

        for k in range(_n - 1):
            for i in range(k + 1, _n):
                fator = U[i][k] / U[k][k]   # fator

                L[i][k] = fator

                for j in range(k, _n):
                    U[i][j] = U[i][j] - fator * U[k][j]  # eliminando [A]

        return L, U
  
    L, U = eliminacao(A)

    # --- Substituições ===============================
    # Resolvendo {L} * {D} = {B}  (substituição PROGRESSIVA, de cima para baixo)
    D[0] = B[0] / L[0][0]
    for i in range(0, n):
        soma = B[i]
        
        for j in range(i):                # j vai de 0 até i-1
            soma -= L[i][j] * D[j]
        D[i] = soma / L[i][i]

    # Resolvendo {U} * {X} = {D}  (substituição REGRESSIVA, de baixo para cima)
    X[n - 1] = D[n - 1] / U[n - 1][n - 1]
    for i in range(n - 2, -1, -1):      # vai até i=0 inclusive
        soma = D[i]
        for j in range(i + 1, n):
            soma -= U[i][j] * X[j]
        X[i] = soma / U[i][i]
        
        
    return np.array(X)


def determinando_coeficientes(vec_A, vec_B) -> np.array:
    '''
    Determina os coeficientes de "a" e "b", calculando uma matriz A a partir de um
    vetor de entrada vec_A. Para determinar o coeficiente, após os cálculos para A,
    o sistema usa o solver_LU() para determinar o vetor dos coeficientes.
    '''
    #vec_A = vec_A[:3]   # pegando só os primeiros 3 termos para fazer o fit
    #vec_B = vec_B[:3]

    def gerando_matriz_A(x):
        A = np.zeros((3, 3))

        A[0][0] = len(x)
        A[0][1] = A[1][0] = np.sum(x)

        for i in range(len(x)):
            A[0][2] += x[i] ** 2
            A[1][2] += x[i] ** 3
            A[2][2] += x[i] ** 4

        A[1][1] = A[2][0] = A[0][2]
        A[2][1] = A[1][2]

        return A


    def gerando_vetor_B(x, y):
        B = np.zeros(3)
        for i in range(0, len(y)):
            B[0] += y[i]
            B[1] += x[i] * y[i]
            B[2] += x[i] ** 2 * y[i]
        return np.array(B)

    matriz_normal = gerando_matriz_A(vec_A)
    vetor_normal = gerando_vetor_B(vec_A, vec_B)
    X = solver_LU(matriz_normal, vetor_normal)

    return np.array(X)


# Calculando os vetores dos coeficientes:
a = determinando_coeficientes(Q, Qt)
b = determinando_coeficientes(Q, Pb)

print(f'a0 = {a[0]}     a1 = {a[1]}     a2 = {a[2]}\n')
print(f'b0 = {b[0]}     b1 = {b[1]}     b2 = {b[2]}\n')



"""
Questão 2 - Otimização
Implementação: Razão Áurea
"""

def func(coef, x):
    """
    coef é o vetor dos coeficientes da função de segunda ordem
    x é o input da função
    Retorna o valor de y da tal função
    """

    return coef[0] + coef[1] * x + coef[2] * x ** 2


def L(x):
    # Retorna o valor da função lucro
    return 0.6 * func(a, x) - 3.0 * func(b, x)


# Resolvendo por Razão Áurea
def razao_aurea():
    R = (-1.0 + 5.0 ** 0.5) / 2.0  # Razão áurea

    xu = 4.0      # valor superior do intervalo de chute
    xl = 1.0      # valor inferior do intervalo de chute

    x = 0
    while xu - xl >= err:
        l1 = (xu - xl) * R
        x2 = xu - l1
        x1 = xl + l1

        if L(x2) > L(x1):
            xu = x1
        else:
            xl = x2
        x += 1

    avg = (xu + xl) / 2
    print(f'Vazão ótima (Q*) = {avg:.5f} m³/h e o lucro máximo é L = $ {L(avg):.2f}\n')
    print(f'Potência térmica útil (Qt) = {func(a, avg):.5f} kW\n')
    print(f'Potência cons. pela bomba (Pb) = {func(b, avg):.5f} kW\n')

razao_aurea()



"""
Questão 3 - Determinando Raízes
Implementação: Newton-Raphson modificado
"""

n = 100     # Nº de iterações desejadas por chute
tol_raiz = 1e-3  # Tolerância para considerar duas raízes iguais
guesses = [-2.0, 0.0, 2.0, 3.0, 5.0, 8.0]  # Chutes iniciais
rootsmod = []   # raízes já encontradas


def func_d(coef, x):
    """
    1ª derivada de func(). Como as duas funções Pt e Qd são polinomiais,
    podemos usar uma func_d comum a elas
    """

    return coef[1] + 2 * coef[2] * x


def func_dd(coef, x):
    """
    2ª derivada de func(). Como as duas funções Pt e Qd são polinomiais,
    podemos usar uma func_d comum a elas
    """

    return 2 * coef[2]


def L_d(x):
    # Retorna o valor da derivada 1ª da função lucro
    return 0.6 * func_d(a, x) - 3.0 * func_d(b, x)


def L_dd(x):
    # Retorna o valor da derivada 2ª da função lucro
    return 0.6 * func_dd(a, x) - 3.0 * func_dd(b, x)


def nr_mod(x):  # Newton-Raphson modificado
    f = L(x)
    fd = L_d(x)
    fdd = L_dd(x)

    denom = fd ** 2 - f * fdd

    if denom == 0:
        denom += 0.1
  
    return x - f * fd / denom


def is_root(y):	 # Verifica se entra na margem de erro
  if err - abs(y) > 0:
    return True
  else:
    return False


def found(x, root_list):    
    # Evita considerar pontos próximos aos que já foram obtidos

  for r in root_list:
    if abs(x - r) < tol_raiz:
      return True
  return False


for guess in guesses:
  xmod = guess

  for i in range(n):
    ymod = L(xmod)

    if is_root(ymod) and not found(xmod, rootsmod):
        """
        Como vazão não pode ser negativa aqui, ficamos com os valores positivos
        """
        if xmod > 0:
            print(f'Vazão crítica (Qc) = {xmod:.3f} m³/h\n')
        rootsmod.append(xmod)
        break
  
    else:
      xmod = nr_mod(xmod)


"""
Questão 4 - Geração do arquivo .dat
"""

from pathlib import Path

def salvar_historico(nome_arquivo, historico):
    """
    Salva o histórico das contas para a tabela em .dat

    As colunas são:

        Q | Qt ajustado | Pb ajustado | Lucro

    :param nome_arquivo: str ou Path
        Nome do arquivo que será criado.

    :param historico: numpy.ndarray
        Matriz retornada pelas contas.

    :return: Path
        Caminho do arquivo criado.
    """

    caminho = Path(nome_arquivo)

    formatos = [
        "%.2e",     # Q
        "%.5e",    # Qt ajustado
        "%.5e",    # Pb ajustado
        "%.5e",    # Lucro
    ]

    np.savetxt(
        caminho,
        historico,
        fmt=formatos,
        header="Q Qt_ajustado Pb_ajustado Lucro"
    )

    return caminho


def calculando_e_salvando(begin, end, steps):
    aQ = begin    # Q atual
    historico = []

    while aQ <= end:
        historico.append(
            [
                aQ,
                func(a, aQ),
                func(b, aQ),
                L(aQ)
            ]
        )
        
        aQ += 0.1

    arquivo_tabela = salvar_historico(
        "resultado.dat",
        historico
    )


calculando_e_salvando(1.0, 4.0, 0.1)

