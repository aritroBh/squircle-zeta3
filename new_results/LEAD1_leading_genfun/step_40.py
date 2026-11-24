import mpmath, json
mpmath.mp.dps = 80
from mpmath import mpf, sqrt, log, quad, pi, gamma, zeta
def kap(i,j):  # ∫ log^i((1+u)/2) log^j((1-u)/2) /[(1-u^2) sqrt(1+u^2)] du
    return quad(lambda u: (log((1+u)/2)**i * log((1-u)/2)**j)/((1-u**2)*sqrt(1+u**2)), [-1,0,1])
def jj(i,j):
    return quad(lambda u: u*(log((1+u)/2)**i * log((1-u)/2)**j)/((1-u**2)*sqrt(1+u**2)), [-1,0,1])
def kapp(i,j):
    return quad(lambda u: (log((1+u)/2)**i * log((1-u)/2)**j)/((1+u)*sqrt(1+u**2)), [-1,0,1])
from math import comb
# sqrt2 sector per d
for d in range(2,9):
    if d%2==1:
        S = sum((-1)**(j+1)*comb(d,j)*jj(d-j,j) for j in range(0,(d-1)//2+1))
    else:
        S = -sum((-1)**(j+1)*comb(d,j)*kap(d-j,j) for j in range(0,d//2)) - mpf('0.5')*(-1)**(d//2+1)*comb(d,d//2)*kapp(d//2,d//2)
    print(f"d={d}: sqrt2-sector = {mpmath.nstr(S,12)}")
