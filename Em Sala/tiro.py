"""
Cálculo Numérico Aplicado - UnB // FT

Atividade do dia 18/06/26:
Método do Tiro
"""
import numpy as np

# Entradas
l = 0.20    # [m] - 20cm
k = 16.2    # [W/m.K] - aço inox, aprox.
qd = 6.37e4 # [W/m³] - 1W pela barra
n = int(input('Número de pontos da malha de simulação:\n> '))

if n <= 1:
  raise(ValueError('O número de pontos deve ser > 1.0 e inteiro.'))

# Condições de contorno
T0 = 25.0  # [ºC] - ambas as extremidades à temp. ambiente
Tl = 25.0  # [ºC]

"""
Sistema versão de primeira ordem:
dT/dx = f
df/dx = -q/k
"""
# Método do tiro
tol = 1e-8  # tolerância

def integrate_f(f_vec, spc, coef) -> np.array:
  """
  f -> vetor vazio em cujos espaços serão alocadas as derivadas calculadas
  apenas o primeiro termo não é nulo, pois é o valor inicial chutado.
  spc -> deltaX entre os pontos da malha (determina os valores de "x")
  """
  num = len(f_vec)
  init_val = f_vec[0]
  for i in range(1, num):
    f_vec[i] = init_val + spc * coef
    init_val = f_vec[i]
    
  return np.array(f_vec)


def integrate_T(vec_T, spc, f_vec) -> np.array:
  """
  vec_T -> o vetor com todos os valores de T para a malha (valores de "y")
  spc -> deltaX entre os pontos da malha (determina os valores de "x")
  f_vec -> vetor das derivadas dT/dx ao longo da barra
  """
  num = len(vec_T)
  init_val = vec_T[0]
  for i in range(1, num):
    vec_T[i] = init_val + spc * f_vec[i - 1]
    init_val = vec_T[i]
    
  return np.array(vec_T)


def euler_integr(f0, T_init):
  f = np.zeros(n)
  f[0] = f0
  T = np.zeros(n)
  T[0] = T_init

  spacing = l / (n - 1)  # deltaX
  coef = - qd / k
  f = integrate_f(f, spacing, coef)  # guarda os valores de dT/dx
  T = integrate_T(T, spacing, f)
  
  return T


f_old = 2000.0
f_new = 6000.0

# Usando f_old:
T = euler_integr(f_old, T0)
err_old = T[-1] - Tl
# print(f'Erro de integração do f_old: {err_old}')

# Usando f_new:
T = euler_integr(f_new, T0)
err_new = T[-1] - Tl
# print(f'Erro de integração do f_new: {err_new}')

# Método da secante
it = 0  # iterações da secante
while abs(err_new) > tol:
  f_next = f_new - err_new * (f_new - f_old) / (err_new - err_old)
  f_old = f_new
  err_old = err_new
  f_new = f_next

  T = euler_integr(f_new, T0)
  err_new = T[-1] - Tl
  it += 1

print(f'''
Valor estimado para f(0) = {f_new}
Temp. obtida em T(L) = {T[-1]}
Erro final: {err_new}
Iterações realizadas: {it}
''')
