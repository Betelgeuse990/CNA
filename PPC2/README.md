# PPC2 - Determinação de raízes de polinômios pelo método de Bairstow

## Resumo operacional

Esta prática tem como objetivo implementar computacionalmente o **método de Bairstow** para a determinação das raízes de polinômios de grau arbitrário, incluindo raízes reais e complexas conjugadas.

A implementação foi utilizada para:

1. determinar raízes de polinômios de teste;
2. validar o algoritmo com um polinômio de sétima ordem construído a partir de raízes conhecidas;
3. analisar a convergência do método para diferentes valores iniciais de $r_0$ e $s_0$;
4. resolver o polinômio característico obtido na APC2;
5. interpretar fisicamente os autovalores encontrados;
6. realizar uma varredura sistemática do plano $(r,s)$;
7. construir mapas de convergência associados ao fractal de Bairstow.

O código principal desta prática está implementado em notebook, no arquivo `main.ipynb`.

## Estrutura do diretório

```text
PPC2/
├── README.md
├── main.ipynb
└── outputs/
```

### Descrição dos arquivos

- **README.md**  
  Documento de explicação da prática, do algoritmo e da organização do diretório.

- **main.ipynb**  
  Notebook principal contendo a implementação do método de Bairstow, a validação do algoritmo, a resolução do polinômio da APC2 e a geração dos gráficos de convergência.

- **outputs/**  
  Diretório destinado aos gráficos e artefatos gerados durante a execução do notebook, como mapas de iterações e regiões de convergência.

## Formulação do problema

O método de Bairstow busca fatorar iterativamente um polinômio em fatores quadráticos da forma:

<p style="text-align:center;">$x^2 - rx - s$</p>

em que $r$ e $s$ são ajustados iterativamente por meio de um processo análogo ao método de Newton-Raphson.

Dado um polinômio

<p style="text-align:center;">$f(x) = a_0 + a_1x + a_2x^2 + \cdots + a_nx^n$</p>

o algoritmo realiza divisões sintéticas sucessivas e calcula correções $\Delta r$ e $\Delta s$ a partir do sistema:

<p style="text-align:center;">$c_1 \Delta r + c_2 \Delta s = -b_0$</p>

<p style="text-align:center;">$c_2 \Delta r + c_3 \Delta s = -b_1$</p>

Após a convergência, obtém-se um fator quadrático e suas raízes associadas, e o processo é repetido sobre o polinômio deflacionado até que todas as raízes sejam encontradas.

## Dicionário de variáveis principais

| Variável | Significado | Unidade / domínio | Tipo |
|---|---|---|---|
| `poly` | coeficientes do polinômio analisado | coeficientes escalares | `list[float]` |
| `r0` | valor inicial de $r$ | adimensional | `float` |
| `s0` | valor inicial de $s$ | adimensional | `float` |
| `r` | valor corrente do parâmetro $r$ | adimensional | `float` |
| `s` | valor corrente do parâmetro $s$ | adimensional | `float` |
| `dr` | correção aplicada a $r$ em uma iteração | adimensional | `float` |
| `ds` | correção aplicada a $s$ em uma iteração | adimensional | `float` |
| `b` | vetor auxiliar da divisão sintética | adimensional | `list[float]` |
| `c` | vetor auxiliar usado no sistema linear de correção | adimensional | `list[float]` |
| `err` | tolerância numérica adotada | adimensional | `float` |
| `max_iter` | número máximo de iterações permitidas | adimensional | `int` |
| `roots` | raízes encontradas para o polinômio | valores complexos | `list[complex]` |
| `poly_validation` | polinômio de validação do item 2 | coeficientes escalares | `list[float]` |
| `poly_apc2` | polinômio característico obtido na APC2 | coeficientes escalares | `list[float]` |
| `autovalores_apc2` | autovalores obtidos do polinômio da APC2 | valores complexos | `list[complex]` |
| `scan_item6` | estrutura com os resultados da varredura do plano $(r,s)$ | mapas numéricos | `dict` |
| `converged_map` | mapa booleano de convergência | lógico | `numpy.ndarray` |
| `iterations_map` | mapa do número de iterações | contagem | `numpy.ndarray` |
| `family_map` | mapa de classificação das regiões de convergência | classe inteira | `numpy.ndarray` |

## Dependências e bibliotecas

A implementação utiliza apenas bibliotecas compatíveis com as diretrizes da disciplina para operações básicas e visualização:

- `cmath`
- `numpy`
- `matplotlib`
- `pandas` (utilizado em algumas tabelas auxiliares)

Não foram utilizados resolvedores prontos para obtenção direta das raízes do problema.

## Entradas e saídas

### Inputs

As entradas são definidas diretamente no notebook, em células específicas, por meio de:

- coeficientes dos polinômios de teste;
- raízes previamente escolhidas para a validação;
- polinômio característico obtido na APC2;
- valores iniciais $r_0$ e $s_0$;
- intervalos e resolução da malha de varredura no plano $(r,s)$;
- tolerância numérica e número máximo de iterações.

### Outputs

A execução do notebook produz:

- raízes encontradas para os polinômios analisados;
- comparação entre raízes exatas e raízes numéricas no teste de validação;
- autovalores do polinômio característico da APC2;
- interpretação física associada ao sistema dinâmico;
- mapa do número de iterações do método;
- mapa das regiões de convergência do método.

Os gráficos gerados devem ser armazenados no diretório `outputs/`.

## Procedimentos de execução

### Opção 1 - Execução em Jupyter Notebook

1. Abra um terminal no diretório do repositório.
2. Navegue até a pasta `PPC2`.
3. Inicie o Jupyter Notebook.
4. Abra o arquivo `main.ipynb`.
5. Execute as células em ordem, do início ao fim.

Exemplo de comandos:

```bash
cd PPC2
jupyter notebook
```

### Opção 2 - Execução no VS Code

1. Abra o repositório no VS Code.
2. Acesse a pasta `PPC2`.
3. Abra o arquivo `main.ipynb`.
4. Execute as células sequencialmente utilizando o kernel Python configurado no ambiente.

## Validação metodológica

A validação do algoritmo foi realizada por meio da construção de um polinômio de sétima ordem com raízes previamente conhecidas, incluindo raízes reais e um par de raízes complexas conjugadas.

A partir dessas raízes conhecidas, o polinômio foi reconstruído computacionalmente e, em seguida, submetido ao método de Bairstow. As raízes numéricas encontradas foram comparadas com as raízes exatas, e os erros absolutos observados permaneceram dentro da tolerância adotada.

Esse procedimento forneceu evidência de que a implementação recupera corretamente as raízes do polinômio de teste, validando o funcionamento do algoritmo para a prática proposta.

## Resultados principais

### 1. Polinômio característico da APC2

O polinômio estudado na APC2 foi:

<p style="text-align:center;">$P(\lambda) = 4.44\lambda^4 + 40.205\lambda^3 + 9863.64\lambda^2 + 37417\lambda + 4537400$</p>

### 2. Autovalores obtidos

Os autovalores numéricos encontrados foram:

<p style="text-align:center;">$\lambda_1 = -3.147960258195 - 39.129231982788j$</p>
<p style="text-align:center;">$\lambda_2 = -3.147960258195 + 39.129231982788j$</p>
<p style="text-align:center;">$\lambda_3 = -1.379629831895 - 25.714949605048j$</p>
<p style="text-align:center;">$\lambda_4 = -1.379629831895 + 25.714949605048j$</p>

### 3. Interpretação física

Como todos os autovalores possuem parte real negativa e parte imaginária não nula, conclui-se que o sistema dinâmico estudado é:

- **estável**;
- **amortecido**;
- **oscilatório em dois modos distintos**.

Além disso, a varredura no plano $(r,s)$ mostrou a existência de duas regiões principais de convergência, cada uma associada a um dos pares de autovalores complexos conjugados do polinômio.

## Bibliografia específica

- CHAPRA, Steven C.; CANALE, Raymond P. **Métodos Numéricos para Engenharia**. 5. ed. McGraw-Hill, 2008.
- GONTIJO, Rafael Gabler. **Notas de aula do curso de Cálculo Numérico Aplicado**. Universidade de Brasília, 2026.
- Documentação oficial da linguagem Python.
- Documentação oficial do NumPy.
- Documentação oficial do Matplotlib.

## Observação final

Esta prática foi desenvolvida com foco em clareza, legibilidade e reprodutibilidade, buscando manter a lógica matemática do método numérico explícita ao longo de toda a implementação.
