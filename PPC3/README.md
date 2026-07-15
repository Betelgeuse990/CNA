# PPC3 — Condução de calor transiente unidimensional

## Consulta rápida

| Pergunta | Resposta |
|---|---|
| Qual problema é resolvido? | Evolução transiente da temperatura em um domínio 1D, com e sem geração interna de calor |
| Quais são os métodos? | Diferenças Finitas, esquema implícito no tempo e algoritmo de Thomas |
| Onde está o algoritmo de Thomas? | [`CNA PPC3.ipynb`](CNA%20PPC3.ipynb), função `thomas_algorithm()` |
| Onde estão os solucionadores térmicos? | `solve_heat_1d_no_generation()` e `solve_heat_1d_with_generation()` |
| Como ocorre a validação? | Comparação com série analítica e recuperação do caso sem geração quando `q_dot → 0` |

## Resumo operacional

O notebook resolve a equação de condução transiente unidimensional. A derivada espacial é aproximada por diferenças centrais e a evolução temporal utiliza um esquema implícito, produzindo um sistema tridiagonal em cada passo. Esse sistema é resolvido por uma implementação explícita do algoritmo de Thomas.

São analisados os casos sem geração interna e com geração volumétrica uniforme. O notebook também produz perfis, mapas espaço-tempo e animações.

## Estrutura

```text
PPC3/
├── README.md
└── CNA PPC3.ipynb
```

## Métodos e rotinas

| Método/etapa | Rotina | Finalidade |
|---|---|---|
| Thomas | `thomas_algorithm()` | resolver o sistema tridiagonal |
| Montagem sem geração | `build_implicit_system_no_generation()` | construir diagonais e vetor independente |
| Solução sem geração | `solve_heat_1d_no_generation()` | avançar a temperatura no tempo |
| Solução analítica | `exact_solution_no_generation()` | fornecer referência para validação |
| Montagem com geração | `build_implicit_system_with_generation()` | incorporar o termo fonte |
| Solução com geração | `solve_heat_1d_with_generation()` | simular o campo térmico com fonte |

## Dicionário de variáveis

| Variável | Significado | Unidade | Tipo |
|---|---|---|---|
| `T` / `T_old` | temperaturas nodais | K ou °C | `numpy.ndarray` |
| `alpha` | difusividade térmica | m²/s | `float` |
| `k` | condutividade térmica | W/(m·K) | `float` |
| `h` | coeficiente de convecção | W/(m²·K) | `float` |
| `L` | comprimento do domínio | m | `float` |
| `dt` | passo temporal | s | `float` |
| `T_inf` | temperatura do ambiente | K ou °C | `float` |
| `q_dot` | geração volumétrica de calor | W/m³ | `float` |
| `lower`, `diag`, `upper` | diagonais do sistema | variável | `numpy.ndarray` |
| `rhs` | vetor independente | variável | `numpy.ndarray` |

## Dependências

- Python 3;
- Jupyter;
- NumPy;
- Matplotlib;
- IPython, instalado juntamente com o Jupyter.

```bash
python -m pip install -r requirements.txt
```

## Entradas e saídas

Os parâmetros físicos e numéricos são definidos nas primeiras células do notebook. As saídas incluem históricos de temperatura, perfis espaciais, mapas de calor, comparações entre modelos e animações incorporadas ao notebook.

## Execução

Na raiz do repositório:

```bash
jupyter notebook "PPC3/CNA PPC3.ipynb"
```

Execute as células sequencialmente, do início ao fim.

## Validação metodológica

1. `thomas_algorithm()` é testado com um sistema tridiagonal de solução conhecida.
2. A solução numérica sem geração é comparada com uma solução analítica em série.
3. O termo de geração é reduzido a zero para verificar se o modelo recupera o caso anterior.
4. O refinamento espacial e temporal permite avaliar a estabilidade dos resultados.

## Hipóteses do modelo

- condução unidimensional;
- propriedades térmicas constantes;
- simetria adiabática em uma extremidade;
- convecção na outra extremidade;
- geração interna uniforme quando ativada.

## Referências

- CHAPRA, Steven C.; CANALE, Raymond P. **Métodos Numéricos para Engenharia**. 5. ed. McGraw-Hill, 2008.
- FINK, J. K. Thermophysical properties of uranium dioxide. *Journal of Nuclear Materials*.
- INTERNATIONAL ATOMIC ENERGY AGENCY. *Thermophysical Properties of Materials for Nuclear Engineering*.
