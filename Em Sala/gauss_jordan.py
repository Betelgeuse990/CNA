# Sistema A X = B
A = [
    [3.0, 4.0, 5.0],
    [6.0, 7.0, 8.0],
    [9.0, 10.0, 1.0]
]

B = [43.0, 73.0, 53.0]

n = len(A)

# Matriz aumentada [A | B]
M = [A[i] + [B[i]] for i in range(n)]

for k in range(n):
    # Pivoteamento parcial
    pivot_row = max(range(k, n), key=lambda i: abs(M[i][k]))

    if abs(M[pivot_row][k]) < 1e-12:
        raise ValueError("O sistema não possui solução única.")

    # Coloca a melhor linha na posição do pivô
    M[k], M[pivot_row] = M[pivot_row], M[k]

    # Normaliza a linha do pivô
    pivot = M[k][k]

    for j in range(n + 1):
        M[k][j] /= pivot

    # Zera a coluna do pivô nas demais linhas
    for i in range(n):
        if i != k:
            fator = M[i][k]

            for j in range(n + 1):
                M[i][j] -= fator * M[k][j]

# A última coluna contém a solução
X = [M[i][-1] for i in range(n)]

print("Matriz reduzida:")
for linha in M:
    print(linha)

print("Solução:", X)
