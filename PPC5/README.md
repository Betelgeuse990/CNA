# PPC #5 — Solução Numérica da Equação de Blasius

Programa desenvolvido para resolver numericamente a equação de Blasius:

$$
f''' + \frac{1}{2} f f'' = 0
$$

sujeita às condições de contorno:

$$
f(0)=0,
\qquad
f'(0)=0,
\qquad
f'(\infty)=1.
$$

A solução é obtida pelo Método do Tiro, utilizando o método da Secante para ajustar o parâmetro desconhecido $f''(0)$. O sistema de equações diferenciais é integrado pelo método de Runge-Kutta de quarta ordem.

---

## Objetivos

O programa realiza as seguintes tarefas:

* transforma a equação de Blasius em um sistema de três EDOs de primeira ordem;
* integra o sistema pelo método de Runge-Kutta de quarta ordem;
* utiliza o Método do Tiro para satisfazer a condição $f'(\eta_{\max}) \approx 1$;
* ajusta o parâmetro de tiro $s=f''(0)$ pelo método da Secante;
* determina numericamente o valor convergido de $f''(0)$;
* salva o perfil de solução contendo $\eta$, $f$, $f'$ e $f''$;
* gera os gráficos dos perfis de similaridade;
* calcula $\eta_{99}$, definido por $f'(\eta_{99})=0.99$;
* calcula o coeficiente local de atrito na parede $C_f$;
* compara os resultados numéricos com valores clássicos da literatura.

---

## Formulação matemática

As variáveis auxiliares utilizadas são:

$$
y_1=f,
\qquad
y_2=f',
\qquad
y_3=f''.
$$

Assim, a equação de Blasius é reescrita como:

$$
\frac{dy_1}{d\eta}=y_2,
$$

$$
\frac{dy_2}{d\eta}=y_3,
$$

$$
\frac{dy_3}{d\eta}
==================

-\frac{1}{2}y_1y_3.
$$

As condições iniciais do problema transformado são:

$$
y_1(0)=0,
\qquad
y_2(0)=0,
\qquad
y_3(0)=s,
$$

em que $s=f''(0)$ é o parâmetro ajustado pelo Método do Tiro.

A função erro utilizada é:

$$
E(s)=f'(\eta_{\max})-1.
$$

O processo iterativo é encerrado quando:

$$
|E(s)|<\varepsilon,
$$

em que $\varepsilon$ é a tolerância definida pelo usuário.

---

## Métodos numéricos empregados

### Método do Tiro

Como a condição $f'(\infty)=1$ não é conhecida em $\eta=0$, o valor inicial $f''(0)$ é tratado como um chute.

Para cada valor de $s$:

1. o sistema é integrado de $\eta=0$ até $\eta=\eta_{\max}$;
2. calcula-se o erro $E(s)=f'(\eta_{\max})-1$;
3. o valor de $s$ é atualizado;
4. o procedimento é repetido até satisfazer a tolerância.

### Método da Secante

A atualização do parâmetro de tiro é feita por:

$$
s_{k+1}
=======

## s_k

E(s_k)
\frac{s_k-s_{k-1}}
{E(s_k)-E(s_{k-1})}.
$$

### Método de Runge-Kutta de quarta ordem

O sistema de três EDOs é integrado simultaneamente por RK4. Para cada passo $h=\Delta\eta$, são calculados os vetores:

$$
K_1=hF(\mathbf{y}_i),
$$

$$
K_2=hF\left(\mathbf{y}_i+\frac{K_1}{2}\right),
$$

$$
K_3=hF\left(\mathbf{y}_i+\frac{K_2}{2}\right),
$$

$$
K_4=hF(\mathbf{y}_i+K_3),
$$

e a atualização é dada por:

$$
\mathbf{y}_{i+1}
================

\mathbf{y}_i+
\frac{K_1+2K_2+2K_3+K_4}{6}.
$$

---

## Estrutura do repositório

```text
PPC5-Blasius/
│
├── Blasius_Friction.py
├── README.md
├── requirements.txt
│
└── resultados/
    ├── perfil_blasius.csv
    ├── perfil_f.png
    ├── perfil_f_linha.png
    ├── perfil_f_duas_linhas.png
    └── resumo_resultados.txt
```

A pasta `resultados/` é criada automaticamente quando o programa é executado.

---

## Requisitos

O programa foi desenvolvido em Python 3 e utiliza as bibliotecas:

```text
numpy
matplotlib
```

Crie um arquivo chamado `requirements.txt` com o seguinte conteúdo:

```text
numpy
matplotlib
```

Instale as dependências com:

```bash
pip install -r requirements.txt
```

---

## Como executar

No terminal, navegue até a pasta do projeto e execute:

```bash
python Blasius_Friction.py
```

Em alguns sistemas, pode ser necessário utilizar:

```bash
python3 Blasius_Friction.py
```

---

## Parâmetros numéricos

Os parâmetros utilizados estão definidos na função `main()` do arquivo `Blasius_Friction.py`.

Exemplo de configuração:

```python
# Chutes iniciais para o Método da Secante
s0 = 0.30
s1 = 0.40

# Parâmetros da integração
delta_eta = 1e-3
eta_max = 10.0
tol = 1e-8
max_tiro = 50
```

| Parâmetro   | Descrição                                  |
| ----------- | ------------------------------------------ |
| `s0`        | Primeiro chute para $f''(0)$               |
| `s1`        | Segundo chute para $f''(0)$                |
| `delta_eta` | Passo de integração $\Delta\eta$           |
| `eta_max`   | Limite superior de integração              |
| `tol`       | Tolerância para o erro do Método do Tiro   |
| `max_tiro`  | Número máximo de atualizações pela Secante |

Os valores `s0 = 0.30` e `s1 = 0.40` são boas escolhas iniciais, pois o valor clássico esperado é próximo de:

$$
f''(0)\approx0.332057.
$$

---

## Dados físicos utilizados no cálculo de atrito

O programa calcula o número de Reynolds local:

$$
Re_x=\frac{U_\infty x}{\nu},
$$

e o coeficiente local de atrito:

$$
C_f=
\frac{2f''(0)}
{\sqrt{Re_x}}.
$$

Os parâmetros físicos estão definidos na função `main()`:

```python
u_inf = 5.0
x = 1.0
nu = 1.5e-5
```

| Variável | Descrição                                    | Unidade |
| -------- | -------------------------------------------- | ------- |
| `u_inf`  | Velocidade do escoamento livre               | m/s     |
| `x`      | Distância medida a partir da borda de ataque | m       |
| `nu`     | Viscosidade cinemática do fluido             | m²/s    |

---

## Arquivos gerados

Após uma execução bem-sucedida, o programa gera os arquivos abaixo dentro da pasta `resultados/`.

### `perfil_blasius.csv`

Arquivo contendo o perfil completo da solução numérica.

| Coluna          | Significado                      |
| --------------- | -------------------------------- |
| `eta`           | Variável de similaridade $\eta$  |
| `f`             | Função de similaridade $f(\eta)$ |
| `f_linha`       | Derivada $f'(\eta)=u/U_\infty$   |
| `f_duas_linhas` | Derivada segunda $f''(\eta)$     |

### `perfil_f.png`

Gráfico de:

$$
f(\eta).
$$

### `perfil_f_linha.png`

Gráfico do perfil de velocidade adimensional:

$$
f'(\eta)=\frac{u}{U_\infty}.
$$

Esse gráfico também indica:

* a linha correspondente a $f'(\eta)=1$;
* a linha correspondente a $f'(\eta)=0.99$;
* o valor numérico calculado de $\eta_{99}$.

### `perfil_f_duas_linhas.png`

Gráfico de:

$$
f''(\eta).
$$

### `resumo_resultados.txt`

Arquivo textual contendo:

* valor convergido de $f''(0)$;
* comparação com o valor clássico;
* número de iterações da Secante;
* erro final do Método do Tiro;
* valor de $f'(\eta_{\max})$;
* valor de $\eta_{99}$;
* valor de $C_\delta=\eta_{99}$;
* comparação com a correlação clássica;
* número de Reynolds local;
* coeficiente local de atrito $C_f$.

---

## Resultados esperados

Para uma configuração com:

```python
s0 = 0.30
s1 = 0.40
delta_eta = 1e-3
eta_max = 10.0
tol = 1e-8
```

espera-se obter valores próximos de:

```text
f''(0) ≈ 0.332057
f'(eta_max) ≈ 1.0
eta_99 ≈ 4.9
C_delta ≈ 4.9
```

O valor clássico utilizado para comparação é:

$$
f''(0)\approx0.332057.
$$

Para a espessura de camada limite definida por $u/U_\infty=0.99$, a correlação clássica é:

$$
\frac{\delta}{x}
================

\frac{4.92}{\sqrt{Re_x}}.
$$

Como:

$$
\frac{\delta}{x}
================

\frac{\eta_{99}}{\sqrt{Re_x}},
$$

tem-se:

$$
C_\delta=\eta_{99}.
$$

---

## Possíveis fontes de erro numérico

As principais fontes de erro são:

* passo de integração $\Delta\eta$ muito grande;
* valor de $\eta_{\max}$ insuficiente para que $f'(\eta)$ se aproxime adequadamente de 1;
* tolerância pouco rigorosa no Método do Tiro;
* número máximo de iterações insuficiente;
* chutes iniciais inadequados para o Método da Secante;
* erro associado à interpolação linear usada para obter $\eta_{99}$.

Uma forma de analisar a convergência é repetir a simulação com valores menores de `delta_eta`, por exemplo:

```python
delta_eta = 1e-2
delta_eta = 5e-3
delta_eta = 1e-3
```

e comparar os valores obtidos para:

$$
f''(0),
\qquad
\eta_{99},
\qquad
C_f.
$$

---

## Autor

**Estevão André**

Disciplina: Cálculo Numérico Aplicado
Professor: Prof. Dr. Rafael Gabler Gontijo
