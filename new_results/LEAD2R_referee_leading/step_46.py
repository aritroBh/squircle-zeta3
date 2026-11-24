import sys; sys.path.insert(0,'.')
import mpmath as mp
mp.mp.dps = 50
# A_{p,q} = int_{-1}^1 L1^p L2^q/(1-u) du ; u=2s-1: = int_0^1 (log s)^p (log(1-s))^q/(1-s) ds
# s = e^-r: = (-1)^p int_0^inf r^p (log(1-e^-r))^q/(e^r-1) dr   [eq (8)]
def A_pq_M(p, q):
    return mp.quad(lambda s: mp.log(s)**p*mp.log(1-s)**q/(1-s), [0, mp.mpf('0.5'), 1])
def A_pq_BE(p, q):
    return (-1)**p*mp.quad(lambda r: r**p*mp.log(1-mp.exp(-r))**q/(mp.exp(r)-1), [0, 1, mp.inf])
worst = mp.mpf(0)
for (p,q) in [(1,0),(2,0),(1,1),(2,1),(3,2),(5,3),(4,4),(6,1),(1,5)]:
    r = abs(A_pq_M(p,q)-A_pq_BE(p,q)); worst = max(worst,r)
print("A_{p,q}: M-form == Bose-Einstein form, worst:", mp.nstr(worst,3))
# sanity A_{p,0} = (-1)^p p! zeta(p+1)
from math import factorial
worst2 = max(abs(A_pq_M(p,0) - (-1)**p*factorial(p)*mp.zeta(p+1)) for p in range(1,6))
print("A_{p,0} == (-1)^p p! zeta(p+1), worst:", mp.nstr(worst2,3))
# Lemma 5.2: H(x,y) = int e^{xL1+yL2}/(1-u) du == B(1+x,y); incl coordinator FAIL point
worst3 = mp.mpf(0)
for (x,y) in [(mp.mpf('0.3'),mp.mpf('0.4')), (mp.mpf('0.3'),mp.mpf('1.4')), (mp.mpf('-0.7'),mp.mpf('0.25')), (mp.mpf('0.9'),mp.mpf('0.9')), (mp.mpc('0.2','0.3'), mp.mpc('0.5','-0.2'))]:
    H = mp.quad(lambda s: s**x*(1-s)**y/(1-s), [0, mp.mpf('0.5'), 1])
    B = mp.beta(1+x, y)
    r = abs(H-B); worst3 = max(worst3,r)
    print(f"  H({mp.nstr(x,4)},{mp.nstr(y,4)}) - B = {mp.nstr(r,3)}")
print("worst H==B:", mp.nstr(worst3,3))
with open('progress.log','a') as f:
    f.write("step_14: A_{p,q} both forms agree (worst "+mp.nstr(worst,3)+"); A_{p,0}=(-1)^p p!zeta(p+1) worst "+mp.nstr(worst2,3)+"; H(x,y)=B(1+x,y) incl (0.3,0.4): worst "+mp.nstr(worst3,3)+" dps50\n")
