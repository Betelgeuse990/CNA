# CNA PPC3 — Condução de Calor Transiente 1D

Este repositório contém a solução do **Programa Para Casa #3 (PPC3)** da disciplina de **Cálculo Numérico Aplicado**.
(Desenvolvido por Estevão A. B. Silva)
O projeto implementa uma simulação numérica de condução de calor transiente unidimensional usando:

- método das diferenças finitas;
- esquema implícito no tempo;
- algoritmo de Thomas para sistemas tridiagonais;
- validação com solução analítica no caso sem geração interna;
- simulação com geração interna de calor;
- visualizações gráficas e animações.

## Arquivo principal

O notebook principal é:

```text
CNA PPC3.ipynb
```

Nele estão a explicação teórica, a implementação dos métodos numéricos, as validações e as visualizações.

## Requisitos

Para executar o notebook, é necessário ter Python instalado.

As bibliotecas utilizadas são:

```text
numpy
matplotlib
jupyter
```

## Instalação

Clone o repositório ou baixe os arquivos manualmente.

```bash
git clone <url-do-repositorio>
cd <nome-do-repositorio>
```

Instale as dependências:

```bash
pip install numpy matplotlib jupyter
```

Caso prefira usar um ambiente virtual:

```bash
python -m venv .venv
```

No Windows:

```bash
.venv\Scripts\activate
```

No Linux/macOS:

```bash
source .venv/bin/activate
```

Depois instale as dependências:

```bash
pip install numpy matplotlib jupyter
```

## Como executar

Abra o Jupyter Notebook:

```bash
jupyter notebook
```

Em seguida, abra o arquivo:

```text
CNA PPC3.ipynb
```

Execute as células em ordem, do início ao fim.

## Estrutura do notebook

O notebook está organizado em quatro partes principais.

### 1. Solução sem geração interna

Nesta etapa, é resolvido o problema de condução de calor transiente unidimensional sem geração interna.

A discretização é feita pelo método das diferenças finitas, usando um esquema implícito no tempo. A cada passo de tempo, surge um sistema linear tridiagonal, resolvido pelo algoritmo de Thomas.

Também é realizada uma comparação com a solução analítica para validar a implementação numérica.

### 2. Solução com geração interna

Nesta etapa, é acrescentado o termo de geração volumétrica de calor ao modelo.

O objetivo é verificar como a presença de uma fonte interna de energia altera a distribuição de temperatura ao longo do domínio.

Também é feita uma validação no limite em que a geração interna tende a zero, mostrando que o modelo com geração recupera o comportamento do caso sem geração.

### 3. Representações visuais

Nesta etapa, são testadas diferentes formas de representação visual da solução:

- perfis de temperatura ao longo da posição;
- mapas de calor espaço-tempo;
- comparação entre os casos com e sem geração interna;
- animação da faixa térmica do corpo.

Essas visualizações ajudam a interpretar o comportamento transiente da solução e o efeito físico da geração interna.

### 4. Valores físicos realistas

Na etapa final, são adotados valores físicos simplificados e plausíveis para representar um combustível sólido de reator nuclear, especialmente dióxido de urânio.

O objetivo não é construir uma simulação nuclear completa, mas aplicar o modelo numérico desenvolvido a um cenário fisicamente inspirado.

## Resultados esperados

Ao executar o notebook, são gerados:

- solução numérica para o caso sem geração interna;
- comparação com a solução analítica;
- solução numérica para o caso com geração interna;
- validação do limite de geração interna nula;
- gráficos de perfis de temperatura;
- mapas de calor;
- animações da evolução térmica;
- simulação final com valores físicos simplificados.

## Observações sobre o modelo

O modelo utilizado é unidimensional e simplificado.

As principais hipóteses são:

- condução de calor em uma direção espacial;
- simetria adiabática em uma extremidade;
- convecção na outra extremidade;
- propriedades térmicas constantes;
- geração interna uniforme no caso com termo fonte.

Essas simplificações tornam o problema adequado para o estudo dos métodos numéricos envolvidos.

## Algoritmo de Thomas

O algoritmo de Thomas foi implementado manualmente para resolver sistemas tridiagonais.

Não foram utilizadas bibliotecas prontas para resolver diretamente o sistema linear principal, pois o objetivo do trabalho é aplicar os conceitos numéricos estudados na disciplina.

## Referências

- International Atomic Energy Agency, *Thermophysical Properties of Materials for Nuclear Engineering*.
- International Atomic Energy Agency, *Thermophysical Properties Database of Materials for Light Water Reactors and Heavy Water Reactors*.
- Fink, J. K., *Thermophysical Properties of Uranium Dioxide*, Journal of Nuclear Materials.
