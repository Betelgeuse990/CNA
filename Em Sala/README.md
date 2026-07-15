# Códigos desenvolvidos em sala e laboratório

Este diretório reúne implementações didáticas, comparações entre métodos e códigos utilizados em avaliações. Ele funciona como um catálogo complementar aos PPCs: localize abaixo a situação-problema, abra o arquivo indicado e procure pela função principal informada.

> Os arquivos possuem níveis diferentes de maturidade. Consulte a coluna **Estado** antes de reutilizar uma implementação.

## Índice rápido

| Situação-problema | Método | Arquivo | Onde está a implementação | Estado |
|---|---|---|---|---|
| Resolver um sistema linear geral | Decomposição LU | [`matrizes_L_U.py`](matrizes_L_U.py) | `eliminacao()` e blocos de substituição | Executável |
| Resolver um sistema pela redução da matriz aumentada | Gauss-Jordan | [`gauss_jordan.py`](gauss_jordan.py) | laços de eliminação e substituição | Executável |
| Encontrar raízes simples e múltiplas | Newton-Raphson comum e modificado | [`newton_raphson.py`](newton_raphson.py) | `nr()` e `nr_mod()` | Executável |
| Comparar Müller e Secante sem derivadas | Müller e Secante | [`muller_vs_secante.py`](muller_vs_secante.py) | `muller()` e `secante()` | Executável |
| Maximizar uma função escalar unimodal | Razão Áurea e Interpolação Quadrática | [`otimizacao_ra_vs_iquad.py`](otimizacao_ra_vs_iquad.py) | `razao_aurea()` e `inter_quad()` | Executável |
| Resolver um PVC de condução 1D | Tiro + Secante + Euler | [`tiro.py`](tiro.py) | `euler_integr()` e laço da Secante | Executável |
| Resolver a equação de Blasius | Tiro + Secante + RK4 | [`Blasius.py`](Blasius.py) | `step_rk4()`, `integrate_rk4()` e laço da Secante | Versão didática |
| Resolver Laplace 2D em uma aleta | Diferenças Finitas + Liebmann | [`liebmann.py`](liebmann.py) | `liebmann()` | Em desenvolvimento |
| Ajustar dados, otimizar lucro e encontrar vazões críticas | Mínimos quadrados + LU + Razão Áurea + Newton modificado | [`p1.py`](p1.py) | `determinando_coeficientes()`, `solver_LU()`, `razao_aurea()` e `nr_mod()` | Avaliação P1 |

## Descrição dos códigos

### `matrizes_L_U.py` — decomposição LU

- **Problema:** resolver `A @ X = B` para uma matriz quadrada geral.
- **Método:** fatoração `A = L @ U`, seguida de substituição progressiva em `L @ D = B` e regressiva em `U @ X = D`.
- **Entrada:** tamanho do sistema definido em `size`; matriz e vetor são gerados no próprio script.
- **Saída:** sistema original, vetor-solução e erro médio da reconstrução `L @ U`.
- **Execução:** `python matrizes_L_U.py`.

### `gauss_jordan.py` — eliminação de Gauss-Jordan

- **Problema pretendido:** resolver um sistema linear por operações elementares.
- **Método:** eliminação sobre a matriz e substituição regressiva.
- **Estado:** o arquivo atualmente gera `IndexError`, pois usa o número de elementos da matriz como ordem do sistema e não inicializa corretamente o vetor-solução. Não deve ser usado como referência pronta antes da correção.
- **Execução para diagnóstico:** `python gauss_jordan.py`.

### `newton_raphson.py` — raízes simples e múltiplas

- **Problema:** localizar as raízes do polinômio `x³ - 6x² + 9x - 4` a partir de diferentes chutes.
- **Métodos:** Newton-Raphson em `nr()` e Newton-Raphson modificado em `nr_mod()`.
- **Objetivo do teste:** comparar a quantidade de iterações, especialmente para a raiz múltipla.
- **Saída:** raízes, chutes e número de iterações no console.
- **Execução:** `python newton_raphson.py`.

O arquivo [`newton_raphson_mod.py`](newton_raphson_mod.py) é atualmente uma cópia idêntica de `newton_raphson.py`. Para consulta, prefira `newton_raphson.py`.

### `muller_vs_secante.py` — comparação de métodos sem derivadas

- **Problema:** encontrar uma raiz de `x³ - 13x - 12` sem utilizar derivadas.
- **Métodos:** Müller em `muller()` e Secante em `secante()`.
- **Objetivo do teste:** comparar convergência e número de iterações a partir de um mesmo chute.
- **Saída:** raiz aproximada, iterações e resíduo no console.
- **Execução:** `python muller_vs_secante.py`.

### `otimizacao_ra_vs_iquad.py` — otimização escalar

- **Problema:** determinar a resistência de carga que maximiza a potência de um circuito.
- **Métodos:** Razão Áurea em `razao_aurea()` e Interpolação Quadrática em `inter_quad()`.
- **Entrada:** parâmetros elétricos e intervalo de busca definidos no início do arquivo.
- **Saída:** resistência ótima, potência máxima e número de iterações.
- **Execução:** `python otimizacao_ra_vs_iquad.py`.

### `tiro.py` — problema de valor de contorno

- **Problema:** calcular a distribuição de temperatura em uma barra com geração volumétrica e temperaturas prescritas nas duas extremidades.
- **Métodos:** transformação em sistema de primeira ordem, integração por Euler, Método do Tiro e atualização do chute pela Secante.
- **Rotinas:** `integrate_f()`, `integrate_T()` e `euler_integr()`.
- **Entrada:** número de pontos solicitado no terminal; propriedades e condições de contorno definidas no script.
- **Saída:** derivada inicial estimada, temperatura final, erro e iterações.
- **Execução:** `python tiro.py`.

### `Blasius.py` — versão didática do solucionador de Blasius

- **Problema:** satisfazer as condições de contorno da equação de Blasius.
- **Métodos:** sistema de três EDOs, RK4, Método do Tiro e Secante.
- **Rotinas:** `sistema_blasius()`, `coefs_rk4()`, `step_rk4()` e `integrate_rk4()`.
- **Saída:** estado final, chute convergido para `f''(0)` e erro do contorno.
- **Execução:** `python Blasius.py`.
- **Versão completa:** consulte [`../PPC5/Blasius_Friction.py`](../PPC5/Blasius_Friction.py).

### `liebmann.py` — condução bidimensional em uma aleta

- **Problema:** resolver a equação de Laplace em regime permanente com temperatura prescrita nas quatro fronteiras.
- **Método:** Diferenças Finitas e Liebmann, equivalente ao Gauss-Seidel aplicado à malha.
- **Rotinas:** `liebmann()`, `mesh_study()` e `post_process()`.
- **Saídas:** dados `x, y, T`, mapa térmico, isotermas, linha central e estudo de malha.
- **Execução:** `python liebmann.py`.
- **Observação:** malhas refinadas podem atingir `max_it` ou exigir tempo elevado; confira a mensagem de convergência antes de usar os resultados.

### `p1.py` — solução integrada da primeira avaliação

- **Problema:** ajustar dados de desempenho, determinar a vazão de maior lucro e localizar vazões críticas.
- **Métodos:** equações normais de mínimos quadrados, decomposição LU, Razão Áurea e Newton-Raphson modificado.
- **Rotinas:** `solver_LU()`, `determinando_coeficientes()`, `razao_aurea()` e `nr_mod()`.
- **Saída:** coeficientes ajustados, vazão ótima, lucro, vazões críticas e `resultado.dat`.
- **Execução:** `python p1.py`.

## Dependências

- Python 3;
- NumPy para `matrizes_L_U.py`, `tiro.py`, `Blasius.py`, `liebmann.py` e `p1.py`;
- Matplotlib para o pós-processamento de `liebmann.py`.

Instale as dependências comuns a partir da raiz:

```bash
python -m pip install -r ../requirements.txt
```
