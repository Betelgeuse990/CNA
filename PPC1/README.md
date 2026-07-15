# PPC1 — Sedimentação de uma esfera com RK4

## Consulta rápida

| Pergunta | Resposta |
|---|---|
| Qual problema é resolvido? | Evolução da velocidade adimensional de uma esfera em sedimentação, com e sem pequena contribuição inercial |
| Qual é o método principal? | Runge-Kutta de quarta ordem (RK4) |
| Onde está o método? | [`main.py`](main.py), funções `ks()`, `rk4()` e `simulate()` |
| Como a solução é validada? | Comparação com soluções analíticas para `Re → 0` e `Re ≠ 0` |
| O que é gerado? | Cinco gráficos na pasta [`results/`](results/) e erros máximos no console |

## Resumo operacional

O programa integra a EDO adimensional da velocidade de uma esfera em regime de baixo Reynolds:

$$
\frac{dv}{dt}=\frac{1}{St}\left(1-v-\frac{3}{8}Re_s v^2\right).
$$

O RK4 é implementado explicitamente e aplicado a diferentes números de Stokes, passos temporais e números de Reynolds da partícula. A prática investiga precisão, refinamento temporal e efeito da parcela inercial.

## Estrutura

```text
PPC1/
├── README.md
├── main.py
└── results/
    ├── item1_st_combined.png
    ├── item2_h_combined.png
    ├── item3_re_nonzero.png
    ├── item4_validation.png
    └── item5_res_sweep.png
```

## Dicionário de variáveis

| Variável | Significado | Unidade/domínio | Tipo |
|---|---|---|---|
| `st` | número de Stokes | adimensional | `float` |
| `res` | número de Reynolds da partícula | adimensional | `float` |
| `v` | velocidade adimensional | adimensional | `float` |
| `t` | tempo adimensional | adimensional | `float` |
| `h` ou `passo` | passo temporal | adimensional | `float` |
| `tf` | tempo final | adimensional | `float` |
| `results` | pares ordenados `(t, v)` | adimensional | `list[tuple]` |

## Dependências

- Python 3;
- Matplotlib.

Instalação pela raiz do repositório:

```bash
python -m pip install -r requirements.txt
```

## Entradas e saídas

As entradas são definidas nas funções `item1_compare_re0()` a `item5_res_sweep()`: valores de `St`, `Re_s`, passo temporal e intervalo de integração.

O programa imprime os erros máximos e salva os gráficos na pasta `results/`. A pasta é criada automaticamente por `save_plot()`.

## Execução

Execute dentro do diretório para que a pasta `results/` seja criada no local correto:

```bash
cd PPC1
python main.py
```

## Validação metodológica

- Para `Re → 0`, `analytic_re0()` fornece a solução analítica utilizada como referência.
- Para o caso inercial, `analytic_inertia()` fornece a segunda referência.
- `max_error()` calcula o maior erro absoluto entre os valores numéricos e analíticos.
- O refinamento temporal demonstra a redução do erro quando `h` diminui.

## Referências

- SOBRAL, Y. D.; OLIVEIRA, T. F.; CUNHA, F. R. *On the unsteady forces during the motion of a sedimenting particle*. Powder Technology, v. 178, p. 129-141, 2007.
- CHAPRA, Steven C.; CANALE, Raymond P. **Métodos Numéricos para Engenharia**. 5. ed. McGraw-Hill, 2008.
