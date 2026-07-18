# PPC #6 — Condução Bidimensional de Calor em uma Aleta

## Consulta rápida

| Pergunta | Resposta |
|---|---|
| Qual problema é resolvido? | Condução bidimensional de calor em regime permanente em uma aleta retangular |
| Qual método de discretização é usado? | Método das Diferenças Finitas |
| Quais métodos resolvem o problema? | Eliminação de Gauss, Liebmann e Liebmann com relaxação |
| Onde está o código? | [`Liebmann-Gauss-Seidel.py`](Liebmann-Gauss-Seidel.py) |
| Como ocorre a validação? | Comparação entre os três métodos e com a solução analítica unidimensional |
| O que é gerado? | Arquivos `.dat` e gráficos armazenados em [`resultados/`](resultados/) |

## Resumo operacional

O programa resolve a equação de Laplace para uma aleta retangular:

$$
\frac{\partial^2 T}{\partial x^2}
+
\frac{\partial^2 T}{\partial y^2}
=0,
$$

considerando:

- temperatura prescrita na base, em $x=0$;
- convecção nas superfícies $y=0$, $y=H$ e na extremidade $x=L$;
- regime permanente;
- propriedades constantes;
- ausência de geração interna de calor.

A discretização por diferenças finitas gera um sistema linear, resolvido por três métodos:

1. Eliminação de Gauss com pivoteamento parcial;
2. Método de Liebmann, equivalente ao Gauss–Seidel;
3. Método de Liebmann com sobre-relaxação.

## Formulação numérica

Para os nós internos:

$$
\frac{T_{i+1,j}-2T_{i,j}+T_{i-1,j}}{\Delta x^2}
+
\frac{T_{i,j+1}-2T_{i,j}+T_{i,j-1}}{\Delta y^2}
=0.
$$

Nas superfícies expostas, utiliza-se a condição de Robin:

$$
-k\frac{\partial T}{\partial n} = h(T-T_\infty)
$$

No método relaxado, a temperatura é atualizada por:

$$
T_{i,j}^{novo} = \omega T_{i,j}^{calculado} + (1-\omega)T_{i,j}^{antigo}
$$

A convergência ocorre quando a maior alteração de temperatura em uma iteração satisfaz:

$$
\max\left|T^{novo}-T^{antigo}\right|
\leq \varepsilon.
$$

## Solução analítica

A temperatura ao longo da linha central é comparada com a solução clássica unidimensional para uma aleta com extremidade convectiva:

$$
\theta(x)=
\frac{
\cosh[m(L-x)]
+
\frac{h}{mk}\sinh[m(L-x)]
}{
\cosh(mL)
+
\frac{h}{mk}\sinh(mL)
},
$$

em que:

$$
\theta=
\frac{T-T_\infty}{T_b-T_\infty},
\qquad
m=
\sqrt{\frac{hP}{kA_c}}.
$$

Foi considerada largura unitária, com $A_c=H$ e $P=2$, representando a troca de calor pelas superfícies superior e inferior.

## Dicionário de variáveis principais

| Variável | Significado | Unidade |
|---|---|---|
| `L` | comprimento da aleta | m |
| `H` | espessura da aleta | m |
| `k` | condutividade térmica | W/(m·K) |
| `h` | coeficiente convectivo | W/(m²·K) |
| `Tb` | temperatura da base | °C |
| `Tamb` | temperatura ambiente | °C |
| `nodes_x`, `nodes_y` | número de nós da malha | adimensional |
| `tol` | tolerância iterativa | °C |
| `omega` | fator de relaxação | adimensional |
| `grid` | campo de temperaturas | °C |

## Dependências

O programa utiliza:

```text
numpy
matplotlib
```
