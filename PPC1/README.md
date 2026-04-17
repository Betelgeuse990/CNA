# PPC1 - Sedimentação de uma esfera com RK4

## Resumo

Este trabalho implementa o método de Runge-Kutta de quarta ordem (RK4) para resolver numericamente a equação diferencial que descreve a sedimentação de uma esfera em regime de baixo Reynolds.

O programa foi utilizado para:

1. comparar a solução numérica com a solução analítica no caso assintótico $Re \to 0$ para diferentes valores de $St$;
2. analisar o efeito do passo de tempo $h$ na qualidade da solução;
3. resolver o caso com pequeno efeito inercial no fluido $(Re \neq 0)$;
4. validar a solução numérica com base na solução analítica disponível;
5. comparar o comportamento da solução para diferentes valores de $Re_s$.

## Estrutura

```text
PPC1/
├── README.md
├── main.py
└── results/
```

- `main.py`: código principal da prática.
- `results/`: gráficos gerados durante a execução.

## Dependências

- `numpy`
- `matplotlib`

## Entradas e saídas

### Entradas
As entradas são definidas diretamente no código, incluindo:

- número de Stokes $St$;
- número de Reynolds da partícula $Re_s$;
- passo de tempo $h$;
- tempo final de simulação;
- condição inicial.

### Saídas
O programa gera:

- solução numérica ao longo do tempo;
- comparação com soluções analíticas;
- gráficos para análise dos diferentes itens da prática.

## Execução

No terminal, dentro da pasta `PPC1`, execute:

```bash
python main.py
```

## Validação

A validação foi feita comparando os resultados numéricos com soluções analíticas conhecidas em dois casos:

- regime assintótico $Re \to 0$;
- caso com pequeno efeito inercial $(Re \neq 0)$.

A concordância entre as curvas numéricas e analíticas indica que a implementação do RK4 está consistente para o problema estudado.

## Resultados

Os resultados mostram que:

- o método RK4 reproduz corretamente a solução analítica para $Re \to 0$;
- a escolha do passo temporal influencia diretamente a qualidade da solução;
- a inclusão do efeito inercial altera a dinâmica da partícula em relação ao caso assintótico;
- diferentes valores de $Re_s$ produzem desvios visíveis em relação ao regime de Stokes puro.

## Referências

- SOBRAL, Y. D.; OLIVEIRA, T. F.; CUNHA, F. R. *On the unsteady forces during the motion of a sedimenting particles*. Powder Technology, 178 (2007), 129-141.
- CHAPRA, S. C.; CANALE, R. P. *Métodos Numéricos para Engenharia*. 5. ed. McGraw-Hill, 2008.
