import random as r
import numpy as np

size = 3    # Tamanho da matriz desejada

err = 1e-10  # Erro tolerável para considerar algo como 0
err_eliminacao = 0.0  # Erro da eliminação progressiva

# --- Matriz 10x10 ================================
A = np.array([
    [r.random() for i in range(size)] for j in range(size)  # Matriz A
])

#A = np.array([  # Prova real se o X = [1, 2, 3] será encontrado
#  [1, 2, 5],
#  [2, 2, 2],
#  [3, 10, 1]
#])

B = np.array([r.random() for i in range(size)])  # Vetor B
# B = np.array([20, 12, 26])  # corrigido: era 22, mas para X=[1,2,3] dá 20

def initial_sys():
  id = 0
  for line in A:
    for idx, el in enumerate(line):
      print(f'{el}*x{idx + 1}', end='')
      if idx + 1 != len(line):
        print(' + ', end='')
      else:
        print(f' = {B[id]}')
    id += 1

D = np.zeros(size)  # Vector D
X = np.zeros(size)  # Vector X

n = len(A)

initial_sys()  # Imprime o sistema linear

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
_A = L @ U  # Teoricamente, matriz A

# --- Calculando Erro =============================
tot_err = 0.0

for j in range(n):
  for i in range(n):
    tot_err += abs(A[i][j] - _A[i][j])

err_eliminacao = tot_err / n ** 2.0

if err_eliminacao > err:
  raise ValueError(f'\nOps... Erro na operação excede cota máxima definida: erro máx: {err} // erro obtido: {err_eliminacao}.')

# --- Substituições ===============================
# Resolvendo {L} * {D} = {B}  (substituição PROGRESSIVA, de cima para baixo)
D[0] = B[0] / L[0][0]
for i in range(1, n):
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

print('\nResultado:')
for id, line in enumerate(X):  # Resposta final
  print(f'x{id + 1}\t' + str(line))
print('_' * 30)
print(f'\nErro na eliminação:\n{err_eliminacao}')
