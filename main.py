from math import sqrt, exp, pi, erf

def rho_norm(x, mu=0, s=1):
  return 1/sqrt(2*pi*s)*exp(-(x-mu)**2/2/s**2)


def f_norm(x, mu=0, s=1):
  return (1+erf((x-mu)/sqrt(2)/s))/2

def inv_f_norm(p, mu, s, t=0.001):
# сначала стандартизуем
  if mu != 0 or s != 1:
    return mu + s * inv_f_norm(p,0,1,t)
  # ищем в полосе значений -100…100 (можно изменить)
  low_x, low_p = -100.0, 0
  hi_x, hi_p = 100.0, 1
  while hi_x - low_x > t:
    mid_x = (low_x + hi_x)/2
    mid_p = f_norm(mid_x)
    if mid_p < p:
      low_x, low_p = mid_x, mid_p
    elif mid_p > p:
      hi_x, hi_p = mid_x, mid_p
    else:
      break
  return mid_x


n = 100
p = 0.11
nu = n*p
o = sqrt(n*p*(1-p))

left = inv_f_norm(0.025,nu,o)
right = 2*nu-inv_f_norm(0.025,nu,o)

# ошибка первого рода
if round(f_norm(right,nu,o)-f_norm(left,nu,o),2) == 0.95:
#ошибка второго родa
  p1 = 0.95*p
  p2 = 1.05*p
  nu1 = p1*n
  nu2 = p2*n
  o1 = sqrt(n*p1*(1-p1))
  o2 = sqrt(n*p2*(1-p2))
  
  w1 = 1 - round(f_norm(right,nu1,o1)-f_norm(left,nu1,o1),2)
  w2 = 1 - round(f_norm(right,nu2,o2)-f_norm(left,nu2,o2),2)

  print(f"Гипотеза верна по уровню значимости 5% с мощностью проверки {round((w1+w2)*50,1)}%, если стрелок промахнется от {round(left)} до {round(right)} раз.")
else:
  print("Гипотеза не верна")