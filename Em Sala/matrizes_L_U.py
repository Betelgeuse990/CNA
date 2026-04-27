import random as r
import numpy as np

size = 3    # Tamanho da matriz desejada

err = 1e-10  # Erro tolerável para considerar algo como 0
err_eliminacao = 0.0  # Erro da eliminação progressiva

# --- Matriz 10x10 ================================
A = np.array([
    [r.random() for i in range(size)] for j in range(size)  # Matriz A
])

B = np.array([r.random() for i in range(size)])  # Vetor B

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
  U = np.array(A[:])  # A tal matriz "U"
  L = np.array([[1.0 if i == j else 0.0 for j in range(_n)] for i in range(_n)])    # A tal matriz "L", identidade
  
  for k in range(n - 1):
      for i in range(k + 1, n):
        fator = U[i][k] / U[k][k]   # fator
          
        L[i][k] = fator
        
        for j in range(k, n):
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

# --- Substituição regressiva =====================
n -= 1  # Ajustando os índices

# Resolvendo {L} * {D} = {B}
D[n] = B[n] / L[n][n]
for i in range(n - 1, 0, -1):
  soma = B[i]
  for j in range(i + 1, n + 1):
    soma -= L[i][j] * D[j]
  D[i] = soma / L[i][i]

# Resolvendo {U} * {X} = {D}
X[n] = D[n] / U[n][n]
for i in range(n - 1, 0, -1):
  soma = D[i]
  for j in range(i + 1, n + 1):
    soma -= U[i][j] * X[j]
  X[i] = soma / U[i][i]

print('\nResultado:')
for id, line in enumerate(X):  # Resposta final
  print(f'x{id + 1}\t' + str(line))
print('_' * 30)
print(f'\nErro na eliminação:\n{err_eliminacao}')
