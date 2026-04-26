# Algorithm to solve a linear system using the Gauss-Jordan method

A = [
    [3, 4, 5],
    [6, 7, 8],
    [9, 10, 1]
]

n = len(A) ** 2    # Quantidade de elementos na matriz A

X = []

B = [43, 73, 53]

# ----------------- Eliminação progressiva -----------------
for k in range(n - 1):
    for i in range(k, n):
        fator = A[i][k] / A[k][k]   # fator
        for j in range(k, n):
            A[i][j] = A[i][j] - fator * A[k][j]   # eliminando [A]  # ty:ignore[invalid-assignment]
            A[i][k] = 0
        B[i] = B[i] - fator * B[k]               # alterando {B}  # ty:ignore[invalid-assignment]


# ----------------- Substituição regressiva -----------------
X[n] = B[n] / A[n][n]   # começando de baixo
for i in range(n - 1, 0, -1):
    soma = B[i]         # subindo de baixo para cima
    for j in range(i + 1, n + 1):
        soma = soma - A[i][j] * X[j]  # soma
    X[i] = soma / A[i][i]             # substituindo em {X}

