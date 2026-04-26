# Encontrando raízes com Newton-Raphson (e versão modificada, para comparar)
# Aula no lab (2 de Abril, 2026)

err = 1e-7	                          # Erro para considerar raiz
tol_raiz = 1e-3                       # Tolerância para considerar duas raízes iguais
n = 100	                              # Nº de iterações desejadas por chute
guesses = [-2.0, 0.0, 2.0, 3.0, 5.0, 8.0]  # Chutes iniciais

rootsn = []                           # Raízes já encontradas
rootsmod = []

def func(x):	# f(x)
  return x ** 3 - 6 * x ** 2 + 9 * x - 4


def funcd(x):	# f'(x)
  return 3 * x ** 2 - 12 * x + 9


def funcdd(x):	# f''(x)
  return 6 * x - 12


def nr(x):	    # Newton-Raphson normal
  f = func(x)
  fd = funcd(x)

  if fd == 0:
    fd += 1
  
  return x - f / fd


def nr_mod(x):	# Newton-Raphson para raízes múltiplas
  f = func(x)
  fd = funcd(x)
  fdd = funcdd(x)
  denom = fd ** 2 - f * fdd

  if denom == 0:
    denom += 0.1
  
  return x - f * fd / denom


def is_root(y):	 # Verifica se entra na margem de erro
  if err - abs(y) > 0:
    return True
  else:
    return False


def found(x, root_list):
  for r in root_list:
    if abs(x - r) < tol_raiz:
      return True
  return False

print('=-= Newton-Raphson =-=')
for guess in guesses:
  xn = guess
  
  for i in range(n):
    yn = func(xn)
  
    if is_root(yn) and not found(xn, rootsn):
      print(f'Raiz {xn:.3f}  |  It. {i}  |  Chute={guess}')
      rootsn.append(xn)
      break
      
    else:
      xn = nr(xn)
      
    
print('\n\n=-= Newton-Raphson Modificado =-=')
for guess in guesses:
  xmod = guess

  for i in range(n):
    ymod = func(xmod)

    if is_root(ymod) and not found(xmod, rootsmod):
      print(f'Raiz {xmod:.3f}  |  It. {i}  |  Chute={guess}')
      rootsmod.append(xmod)
      break
  
    else:
      xmod = nr_mod(xmod)
