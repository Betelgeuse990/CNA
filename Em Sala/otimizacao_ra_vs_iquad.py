"""
Razão Áurea vs Interpolação Quadrática
Otimizando a potência de um circuito, determinando um valor para uma resistência ótima e a potência max. nesse cenário.
"""

err = 1e-6

V = 80.0  # Tensão na malha
r1 = 8.0  # Ohm
r2 = 12.0
r3 = 10.0

R = (-1.0 + 5.0 ** 0.5) / 2.0  # Razão áurea

val_1, val_2 = 25.0, 0.0

def pwr(Ra) -> float:  # Potência para uma resistência Ra
  if Ra <= 0.0:
    return 0.0
  upper = V * r3 * Ra / (r1 * (Ra + r2 + r3) + r3 * Ra + r3 * r2)
  power = (upper ** 2.0) / Ra
  return power


# Resolvendo por Razão Áurea
def razao_aurea():
  xu = val_1  # intervalo superior de chute
  xl = val_2

  x = 0
  while xu - xl >= err:
    l1 = (xu - xl) * R
    x2 = xu - l1
    x1 = xl + l1

    if pwr(x2) > pwr(x1):
      xu = x1
    else:
      xl = x2
    x += 1

  avg = (xu + xl) / 2
  print(f'Ra_opt = {avg} Ω e P_max = {pwr(avg)} W; it: {x}')


# Resolvendo por interpolação quadrática
def inter_quad():
  x0 = val_2
  x2 = val_1
  x1 = (x0 + x2) / 2

  def next_X(v0, v1, v2):
    f0 = pwr(v0)
    f1 = pwr(v1)
    f2 = pwr(v2)

    upp = (
        f0 * (v1 ** 2 - v2 ** 2)
        + f1 * (v2 ** 2 - v0 ** 2)
        + f2 * (v0 ** 2 - v1 ** 2)
    )

    lwr = (
        2 * f0 * (v1 - v2)
        + 2 * f1 * (v2 - v0)
        + 2 * f2 * (v0 - v1)
    )

    if abs(lwr) < 1e-15:
        return (v0 + v2) / 2

    return upp / lwr
  
  it = 0

  while abs(x2 - x0) >= err:
      x3 = next_X(x0, x1, x2)

      if x3 <= x0 or x3 >= x2:
          x3 = (x0 + x2) / 2

      f1 = pwr(x1)
      f3 = pwr(x3)

      if x3 > x1:
          if f3 > f1:
              # O novo ponto é melhor e está à direita
              x0 = x1
              x1 = x3
          else:
              # O novo ponto é pior e está à direita
              x2 = x3

      else:
          if f3 > f1:
              # O novo ponto é melhor e está à esquerda
              x2 = x1
              x1 = x3
          else:
              # O novo ponto é pior e está à esquerda
              x0 = x3

      it += 1

  print(f'Ra_opt = {x1} Ω e P_max = {pwr(x1)} W; it: {it}')

print('Razão Áurea:')
razao_aurea()
print('Interp. Quadrática:')
inter_quad()

"""
Nos testes, a razão áurea se mostrou bem previsível, enquanto a interpolação quadrática, "temperamental".
A RA demorou por volta de 35 iterações para três intervalos diferentes de valores de val_1 (50.0, 25.0, 20.0), enquanto a interp quad. demorou 88 para o primeiro caso, 21 para o segundo e 13 para o terceiro.
Dessa forma, a quadrática PODE ser muito mais ágil ou muito mais lenta que a da razão áurea com as duas sub as mesmas circunstâncias.
"""
