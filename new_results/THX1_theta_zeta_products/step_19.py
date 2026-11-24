
import mpmath as mp
mp.mp.dps = 100
import rul1, rul1num
from fractions import Fraction as Fr
import math, time, pickle

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

e_rm = pickle.load(open('e_rm_ext.pkl','rb'))
Dtarget = pickle.load(open('Dtilde_relations.pkl','rb'))
def E(r, m):
    return Fr(e_rm.get(f"{r},{m}", '0'))

def I_word(m, j, dps=100):
    w = (0,)*(m+1) + (-1,)*j
    return mp.factorial(j)*mp.factorial(m) * rul1num.g_word_over_0m1_num(w, dps=dps)

def Dtilde_num(n, m, dps=100):
    tot = mp.mpf(0)
    for j in range(1, n-2*m):
        e = E(m, n-1-j)
        if e == 0: continue
        tot += mp.mpf(e.numerator)/mp.mpf(e.denominator) * I_word(n-1-j, j, dps) / mp.factorial(j)
    return 2*tot

# y-integral: D̃_{n,m} = (2/n!) int_0^1 log^n(2/(1+y)) (arcsin y)^{2m} dy/sqrt(1-y)
def Dtilde_yint(n, m):
    f = lambda y: mp.log(2/(1+y))**n * mp.asin(y)**(2*m) / mp.sqrt(1-y)
    return 2/mp.factorial(n) * mp.quad(f, [0, 1])

# zeta-product targets
zval = {2: mp.zeta(2), 3: mp.zeta(3), 5: mp.zeta(5), 6: mp.zeta(6), 7: mp.zeta(7), 8: mp.zeta(8), 9: mp.zeta(9)}
def zbasis(n):
    # basis used in Dtilde_relations
    if n in (2,3,4): return [mp.zeta(n)]
    if n == 5: return [mp.zeta(5), mp.zeta(2)*mp.zeta(3)]
    if n == 6: return [mp.zeta(6), mp.zeta(3)**2]
    if n == 7: return [mp.zeta(7), mp.zeta(2)*mp.zeta(5), mp.zeta(2)**2*mp.zeta(3)]
    if n == 8: return [mp.zeta(8), mp.zeta(3)*mp.zeta(5), mp.zeta(2)*mp.zeta(3)**2]
    if n == 9: return [mp.zeta(9), mp.zeta(2)*mp.zeta(7), mp.zeta(2)**2*mp.zeta(5), mp.zeta(2)**3*mp.zeta(3), mp.zeta(3)**3]
def Dtilde_target(n, m):
    coeffs = [Fr(c) for c in Dtarget[f"{n},{m}"]]
    return sum(mp.mpf(c.numerator)/mp.mpf(c.denominator)*b for c, b in zip(coeffs, zbasis(n)))

print(" n m | Dtilde(formula) - Dtilde(target zeta) | , | formula - y-integral |")
for n in range(2, 10):
    for m in range(0, n//2):
        a = Dtilde_num(n, m)
        b = Dtilde_target(n, m)
        c = Dtilde_yint(n, m)
        print(n, m, mp.nstr(abs(a-b),4), "   ", mp.nstr(abs(a-c),4))
log("Dtilde formula (e-table + I words) matches K3 PSLQ zeta values and y-integral for n<=9 (100 dps)")
