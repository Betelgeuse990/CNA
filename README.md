# Cálculo Numérico Aplicado (CNA)

**Universidade de Brasília (UnB)**  
**Faculdade de Tecnologia (FT)**  
**Departamento de Engenharia Mecânica**  
**Disciplina:** Cálculo Numérico Aplicado  
**Semestre:** 2026/1  
**Discente:** Estevão A. B. Silva

## Sobre o repositório

Este repositório reúne implementações de métodos numéricos desenvolvidas em práticas para casa, exercícios de laboratório e avaliações da disciplina. Sua organização foi pensada como uma biblioteca de consulta: a partir do tipo de problema, o índice abaixo indica o método apropriado, o arquivo que contém a implementação e a principal função a ser consultada.

Os algoritmos numéricos foram implementados explicitamente. NumPy e Matplotlib são utilizados para operações básicas, armazenamento de dados e visualização, sem substituir a lógica dos métodos estudados.

## Índice rápido: problema, método e código

### Sistemas lineares e ajuste de curvas

| Situação-problema | Método | Código | Rotina principal | Estado |
|---|---|---|---|---|
| Resolver um sistema linear geral | Decomposição LU | [`matrizes_L_U.py`](Em%20Sala/matrizes_L_U.py) | `eliminacao()` e substituições triangular progressiva/regressiva | Executável |
| Resolver um sistema linear por redução da matriz aumentada | Gauss-Jordan com pivoteamento parcial | [`gauss_jordan.py`](Em%20Sala/gauss_jordan.py) | normalização da linha pivô e eliminação acima e abaixo da diagonal | Executável |
| Ajustar uma parábola a dados experimentais | Equações normais dos mínimos quadrados + LU | [`p1.py`](Em%20Sala/p1.py) | `determinando_coeficientes()` e `solver_LU()` | Executável |
| Resolver sistemas tridiagonais de condução 1D | Algoritmo de Thomas | [`CNA PPC3.ipynb`](PPC3/CNA%20PPC3.ipynb) | `thomas_algorithm()` | PPC concluído |

### Raízes de equações e polinômios

| Situação-problema | Método | Código | Rotina principal | Estado |
|---|---|---|---|---|
| Encontrar raízes simples de uma função com derivada disponível | Newton-Raphson | [`newton_raphson.py`](Em%20Sala/newton_raphson.py) | `nr()` | Executável |
| Encontrar raízes múltiplas com maior eficiência | Newton-Raphson modificado | [`newton_raphson.py`](Em%20Sala/newton_raphson.py) | `nr_mod()` | Executável |
| Encontrar raízes sem calcular derivadas | Método da Secante | [`muller_vs_secante.py`](Em%20Sala/muller_vs_secante.py) | `secante()` | Executável |
| Comparar aproximação quadrática e secante | Müller e Secante | [`muller_vs_secante.py`](Em%20Sala/muller_vs_secante.py) | `muller()` e `secante()` | Executável |
| Encontrar todas as raízes reais e complexas de um polinômio | Método de Bairstow | [`main.ipynb`](PPC2/main.ipynb) | `bairstow_roots()` | PPC concluído |

### Otimização

| Situação-problema | Método | Código | Rotina principal | Estado |
|---|---|---|---|---|
| Maximizar uma função escalar unimodal em um intervalo | Razão Áurea | [`otimizacao_ra_vs_iquad.py`](Em%20Sala/otimizacao_ra_vs_iquad.py) | `razao_aurea()` | Executável |
| Estimar o máximo por uma parábola interpoladora | Interpolação Quadrática | [`otimizacao_ra_vs_iquad.py`](Em%20Sala/otimizacao_ra_vs_iquad.py) | `inter_quad()` | Executável |
| Maximizar uma função de várias variáveis seguindo o gradiente | Aclive Máximo | [`CNA PPC4.ipynb`](PPC4/CNA%20PPC4.ipynb) | `aclive_maximo()` | PPC concluído |
| Reduzir o zigue-zague da otimização por gradiente | Fletcher-Reeves | [`CNA PPC4.ipynb`](PPC4/CNA%20PPC4.ipynb) | `fletcher_reeves()` | PPC concluído |
| Ajustar curvas e maximizar lucro a partir de dados experimentais | Mínimos quadrados + LU + Razão Áurea | [`p1.py`](Em%20Sala/p1.py) | `determinando_coeficientes()` e `razao_aurea()` | Avaliação P1 |

### Equações diferenciais ordinárias e problemas de contorno

| Situação-problema | Método | Código | Rotina principal | Estado |
|---|---|---|---|---|
| Integrar uma EDO de primeira ordem com alta precisão | Runge-Kutta de quarta ordem | [`main.py`](PPC1/main.py) | `rk4()` e `simulate()` | PPC concluído |
| Simular o movimento gravitacional de dois corpos e comparar conservação de energia | Euler explícito e Leapfrog/Velocity-Verlet | [`euler_vs_leapfrog.py`](Em%20Sala/euler_vs_leapfrog.py) | `integrate_euler()` e `integrate_leapfrog()` | Executável |
| Transformar um problema de contorno em problemas de valor inicial | Método do Tiro + Secante + Euler | [`tiro.py`](Em%20Sala/tiro.py) | `euler_integr()` e laço da Secante | Executável |
| Resolver a equação de Blasius | Tiro + Secante + RK4 | [`Blasius.py`](Em%20Sala/Blasius.py) | `step_rk4()`, `integrate_rk4()` e laço da Secante | Versão didática |
| Resolver Blasius e calcular grandezas da camada limite | Tiro + Secante + RK4 | [`Blasius_Friction.py`](PPC5/Blasius_Friction.py) | `metodo_tiro_secante()` | PPC concluído |

### Salvamento de dados e visualização

| Necessidade | Código de referência | Rotina ou trecho | Saída |
|---|---|---|---|
| Salvar posição, velocidade e energia de uma integração temporal | [`euler_vs_leapfrog.py`](Em%20Sala/euler_vs_leapfrog.py) | `save_results()` | arquivos `.dat` |
| Plotar trajetórias orbitais | [`euler_vs_leapfrog.py`](Em%20Sala/euler_vs_leapfrog.py) | `plot_trajectory()` | PNG |
| Plotar espaço de fase | [`euler_vs_leapfrog.py`](Em%20Sala/euler_vs_leapfrog.py) | `plot_phase_space()` | PNG |
| Comparar conservação de energia entre integradores | [`euler_vs_leapfrog.py`](Em%20Sala/euler_vs_leapfrog.py) | `plot_energy_comparison()` | PNG |
| Salvar tabelas e curvas unidimensionais | [`Blasius_Friction.py`](PPC5/Blasius_Friction.py) | `salvar_perfil_csv()` e `gerar_graficos()` | CSV, TXT e PNG |
| Salvar e recarregar históricos iterativos | [`CNA PPC4.ipynb`](PPC4/CNA%20PPC4.ipynb) | `salvar_historico()` e `plotar_trajetorias()` | DAT e PNG |
| Plotar campos bidimensionais e curvas de nível | [`liebmann.py`](Em%20Sala/liebmann.py) | `post_process()` | DAT e PNG |
| Salvar estudos de relaxação, refinamento de malha e campos 2D | [`Liebmann-Gauss-Seidel.py`](PPC6/Liebmann-Gauss-Seidel.py) | `relaxation_study()`, `mesh_study()` e `post_process()` | DAT e PNG |
| Plotar várias soluções numéricas e analíticas | [`main.py`](PPC1/main.py) | `save_plot()` e funções `item*()` | PNG |

### Equações diferenciais parciais e transferência de calor

| Situação-problema | Método | Código | Rotina principal | Estado |
|---|---|---|---|---|
| Simular condução transiente unidimensional | Diferenças Finitas implícitas + Thomas | [`CNA PPC3.ipynb`](PPC3/CNA%20PPC3.ipynb) | `solve_heat_1d_no_generation()` e `solve_heat_1d_with_generation()` | PPC concluído |
| Resolver a equação de Laplace em uma aleta 2D | Diferenças Finitas + Liebmann (Gauss-Seidel) | [`liebmann.py`](Em%20Sala/liebmann.py) | `liebmann()`, `mesh_study()` e `post_process()` | Executável |
| Resolver a equação de Laplace em uma aleta 2D e comparar métodos de solução | Diferenças Finitas + Gauss + Liebmann + sobre-relaxação | [`Liebmann-Gauss-Seidel.py`](PPC6/Liebmann-Gauss-Seidel.py) | `build_system()`, `gaussian_elimination()`, `liebmann()` e `analytical_solution()` | PPC concluído |

## Práticas para casa

| Diretório | Situação-problema | Métodos principais | Arquivo principal |
|---|---|---|---|
| [`PPC1`](PPC1/) | Sedimentação de uma esfera em baixo Reynolds | RK4 e validação analítica | [`main.py`](PPC1/main.py) |
| [`PPC2`](PPC2/) | Raízes de polinômios e autovalores | Bairstow | [`main.ipynb`](PPC2/main.ipynb) |
| [`PPC3`](PPC3/) | Condução de calor transiente 1D | Diferenças Finitas implícitas e Thomas | [`CNA PPC3.ipynb`](PPC3/CNA%20PPC3.ipynb) |
| [`PPC4`](PPC4/) | Otimização multidimensional | Aclive Máximo, Fletcher-Reeves e Interpolação Quadrática | [`CNA PPC4.ipynb`](PPC4/CNA%20PPC4.ipynb) |
| [`PPC5`](PPC5/) | Equação de Blasius e camada limite | Tiro, Secante e RK4 | [`Blasius_Friction.py`](PPC5/Blasius_Friction.py) |
| [`PPC6`](PPC6/) | Condução de calor permanente 2D em uma aleta | Diferenças Finitas, Gauss, Liebmann e sobre-relaxação | [`Liebmann-Gauss-Seidel.py`](PPC6/Liebmann-Gauss-Seidel.py) |
| [`Em Sala`](Em%20Sala/) | Exercícios, protótipos e avaliações | Vários métodos | [`README.md`](Em%20Sala/README.md) |

## Topologia do repositório

```text
CNA/
├── README.md                  # Índice geral por problema e por método
├── requirements.txt           # Dependências comuns
├── Em Sala/
│   ├── README.md              # Catálogo dos códigos de laboratório
│   └── *.py                   # Implementações didáticas e avaliações
├── PPC1/
│   ├── README.md
│   ├── main.py
│   └── results/
├── PPC2/
│   ├── README.md
│   ├── main.ipynb
│   └── outputs/
├── PPC3/
│   ├── README.md
│   └── CNA PPC3.ipynb
├── PPC4/
│   ├── README.md
│   ├── CNA PPC4.ipynb
│   └── arquivos de resultados
├── PPC5/
│   ├── README.md
│   ├── Blasius_Friction.py
│   ├── requirements.txt
│   └── resultados/
└── PPC6/
    ├── README.md
    ├── Liebmann-Gauss-Seidel.py
    └── resultados/
```

## Preparação do ambiente

Na raiz do repositório, crie um ambiente virtual e instale as dependências comuns:

```bash
python -m venv .venv
```

No Windows:

```bash
.venv\Scripts\activate
```

No Linux ou macOS:

```bash
source .venv/bin/activate
```

Depois:

```bash
python -m pip install -r requirements.txt
```

Entre no diretório da atividade antes de executar scripts que geram arquivos, garantindo que as saídas sejam criadas na pasta correta. Para os notebooks:

```bash
jupyter notebook
```

As instruções específicas, entradas, saídas e validações estão documentadas no README de cada diretório.

## Convenções adotadas

- Cada PPC possui diretório e README próprios.
- Resultados necessários à reprodução permanecem junto à atividade correspondente.
- O diretório `Em Sala` contém exemplos didáticos, protótipos e códigos de avaliação; seu README informa o estado de cada implementação.
- Funções-chave são indicadas nos índices para permitir acesso rápido ao trecho que implementa cada método.
- Bibliotecas prontas não substituem a implementação dos algoritmos numéricos estudados.

## Bibliografia geral

- CHAPRA, Steven C.; CANALE, Raymond P. **Métodos Numéricos para Engenharia**. 5. ed. McGraw-Hill, 2008.
- GONTIJO, Rafael Gabler. **Notas de aula do curso de Cálculo Numérico Aplicado**. Universidade de Brasília, 2026.
- HARRIS, Charles R. et al. Array programming with NumPy. *Nature*, v. 585, p. 357-362, 2020.
- HUNTER, John D. Matplotlib: A 2D graphics environment. *Computing in Science & Engineering*, v. 9, n. 3, p. 90-95, 2007.
