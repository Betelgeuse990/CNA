# Comparação método Müller vs Secante para o polinômio x³ - 13x - 12

erro = 1e-5  # Margem de erro

# Chute inicial
chute_inicial = 6.5  # Testando com um chute maior, como 50, é possível ver a diferença gritante de performance


def func(x) -> float:  # A função em si
  return x ** 3 - 13 * x - 12


def margem(x):
  y = func(x)
  return erro - abs(y)


def muller(x0, err):
  it = 0
  while True:
    x1 = x0 - 0.5
    x2 = x0 - 1.0
    
    h0 = x1 - x0
    h1 = x2 - x1
    
    d0 = (func(x1) - func(x0)) / (x1 - x0)
    d1 = (func(x2) - func(x1)) / (x2 - x1)
  
    a = (d1 - d0) / (h1 + h0)
    b = a * h1 + d1
    c = func(x2)
  
    x0 = x2 - 2 * c / (b + (b ** 2 - 4 * a * c) ** 0.5)
    
    if margem(x0) >= 0:
      break
      
    it += 1
    
  print(f'Müller // Raiz: x={x0:.4f} // It. {it} // Err.: {margem(x0)}')


def secante(x1, err):
  x0 = x1 - 0.5
  
  it = 0
  
  while True:
    x1 = x1 - func(x1) * (x1 - x0) / (func(x1) - func(x0))

    if margem(x1) >= 0:
      break

    it += 1

  print(f'Secante // Raiz: x={x1:.4f} // It. {it} // Err.: {margem(x1)}')
  

muller(chute_inicial, erro)
secante(chute_inicial, erro)
