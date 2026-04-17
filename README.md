# Cálculo Numérico Aplicado (CNA)

Repositório para as atividades da disciplina de Cálculo Numérico Aplicado (CNA).

**Universidade de Brasília (UnB)**  
**Faculdade de Tecnologia (FT)**  
**Departamento de Engenharia Mecânica**  
**Disciplina:** Cálculo Numérico Aplicado  
**Semestre:** 2026/1  
**Discente:** Estevão A.

## Apresentação

Este repositório foi criado para armazenar, organizar e documentar os trabalhos desenvolvidos na disciplina de **Cálculo Numérico Aplicado**. A proposta central é reunir, em um único ambiente versionado, as implementações computacionais, relatórios, gráficos e artefatos complementares produzidos ao longo da disciplina.

O estudo de métodos numéricos é fundamental em engenharia, pois permite resolver computacionalmente problemas matemáticos para os quais nem sempre há solução analítica simples ou viável. Nesse contexto, o uso de um repositório público e bem estruturado favorece a organização, a reprodutibilidade dos resultados e a clareza da lógica implementada em cada prática.

## Estrutura do repositório

A organização do repositório segue uma divisão por prática, de modo que cada atividade possua seu próprio diretório, contendo os códigos, documentações e artefatos associados.

```text
CNA/
├── README.md
├── Em Sala/
├── PPC1/
└── PPC2/
    ├── README.md
    ├── main.ipynb
    └── outputs/
```

### Descrição dos diretórios

- **README.md**  
  Documento principal do repositório. Contém a apresentação geral da disciplina, a topologia do projeto e a bibliografia geral.

- **Em Sala/**  
  Diretório destinado a materiais, códigos e anotações produzidos durante as atividades em sala de aula.

- **PPC1/**  (mexer depois!)
  Diretório da primeira prática/programa para casa. Deve conter o código-fonte, README interno e eventuais arquivos auxiliares de entrada e saída.

- **PPC2/**  
  Diretório da segunda prática/programa para casa. Contém a implementação do método de Bairstow, o notebook principal, gráficos gerados e o README interno específico da prática.

## Tecnologias utilizadas

As implementações presentes neste repositório utilizam, prioritariamente:

- **Python**
- **Jupyter Notebook**
- **NumPy**
- **Matplotlib**
- **Pandas** (quando necessário para organização/tabulação dos resultados)

O uso dessas bibliotecas foi feito com foco em manipulação de dados, operações matemáticas básicas e geração de gráficos, sem recorrer a resolvedores prontos do problema numérico estudado.

## Organização adotada

Cada prática é organizada com foco em:

- clareza do algoritmo implementado;
- documentação interna específica;
- reprodutibilidade dos resultados;
- separação entre código principal e artefatos gerados;
- facilidade de navegação e leitura por terceiros.

Sempre que necessário, gráficos, imagens, relatórios e outros resultados de processamento são armazenados dentro do diretório correspondente à prática.

## Bibliografia geral

- CHAPRA, Steven C.; CANALE, Raymond P. **Métodos Numéricos para Engenharia**. 5. ed. McGraw-Hill, 2008.
- GONTIJO, Rafael Gabler. **Notas de aula do curso de Cálculo Numérico Aplicado**. Universidade de Brasília, 2026.
- HARRIS, Charles R. et al. **Array programming with NumPy**. *Nature*, v. 585, p. 357–362, 2020.
- HUNTER, John D. **Matplotlib: A 2D graphics environment**. *Computing in Science & Engineering*, v. 9, n. 3, p. 90–95, 2007.
  
## Observações finais

Este repositório está em desenvolvimento contínuo ao longo do semestre e será atualizado à medida que novas práticas forem implementadas, refinadas e documentadas.
