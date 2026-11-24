import mpmath, json
from mpmath import mp, mpf, sqrt, pi, log, zeta, polylog, catalan, nstr
mp.dps = 250
sq2 = sqrt(2)
b = sq2-1; a = 1+sq2
L = log(a); lam = log(2)
c8 = mpf('0.5')*sqrt(2-sq2)   # cos(3pi/8) = sin(pi/8)
s8 = mpf('0.5')*sqrt(2+sq2)   # sin(3pi/8)
def Li(n,x): return polylog(n, x)
# L(chi_8,3) via Hurwitz zeta
Lchi = (zeta(3, mpf(1)/8) - zeta(3, mpf(3)/8) - zeta(3, mpf(5)/8) + zeta(3, mpf(7)/8))/512
pts = {'b':b, '-b':-b, 'b2':b**2, '-b2':-b**2, 'half':mpf('0.5'),
       'isq2':1/sq2, 'c1misq2':1-1/sq2, 'tomsq2':2-sq2, 'third':mpf(1)/3}
C = {}
C['pi']=nstr(pi,240); C['sq2']=nstr(sq2,240); C['lam']=nstr(lam,240); C['L']=nstr(L,240)
C['b']=nstr(b,240); C['a']=nstr(a,240); C['catalan']=nstr(catalan,240)
for n in [2,3,4,5,6]: C[f'zeta{n}']=nstr(zeta(n),240)
for k,x in pts.items():
    for n in [2,3,4]:
        C[f'Li{n}_{k}'] = nstr(Li(n,x),240)
C['Lchi8_3'] = nstr(Lchi,240)
C['chi3_b'] = nstr((Li(3,b)-Li(3,-b))/2, 240)
# 135-degree corner constants
C['c8']=nstr(c8,240); C['s8']=nstr(s8,240); C['ln3']=nstr(log(3),240)
for nm,x in [('lnc8',c8),('ln1mc8',1-c8),('ln1pc8',1+c8),('lnc8sq',c8**2),('ln1mc8sq',1-c8**2)]:
    C[nm]=nstr(log(x),240)
for nm,x in [('c8sq',c8**2),('om8',1-c8**2),('r8',(1-c8)/(1+c8)),('b2',b**2),('half',mpf('0.5'))]:
    for n in [2,3]:
        C[f'Li{n}_8_{nm}'] = nstr(Li(n,x),240)
# Lame lambda=2 (golden) and lambda=3 constants
sq5 = sqrt(5); t2 = 2/(1+sq5)
sq10 = sqrt(10); t3 = 3/(1+sq10)
C['t2']=nstr(t2,240); C['t3']=nstr(t3,240)
C['ln5']=nstr(log(5),240); C['ln10']=nstr(log(10),240); C['lnsq5']=nstr(log(sq5),240); C['lnsq10']=nstr(log(sq10),240)
for nm,x in [('t2',t2),('-t2',-t2),('t2sq',t2**2),('t3',t3),('-t3',-t3),('t3sq',t3**2)]:
    for n in [2,3]:
        C[f'Li{n}_{nm}'] = nstr(Li(n,x),240)
json.dump(C, open('ckpt_constants.json','w'))
print("saved", len(C), "constants")
# quick sanity: chi_3 bridge 48 chi3(b) =? 32 sqrt2 Lchi - 3 pi^2 L + 4 L^3
res = 48*(Li(3,b)-Li(3,-b))/2 - (32*sq2*Lchi - 3*pi**2*L + 4*L**3)
print("chi3 bridge residual:", nstr(res, 8))
# sanity: M0(3pi/4) formula
M0 = -(1+c8)*log(1+c8)-(1-c8)*log(1-c8)
print("M0(3pi/4) =", nstr(M0, 40))
print("Li2_b =", nstr(Li(2,b), 50))
print("Lchi =", nstr(Lchi, 50))
with open('progress.log','a') as f:
    f.write("step_29: constants module (250 dps): squircle/135-deg/Lame alphabets saved to ckpt_constants.json; chi3 bridge residual printed\n")
