script = '''import sympy as sp
p = sp.Symbol('p', integer=True, positive=True)
i = sp.I
t = { (1,0): (-1-i)/2, (2,0): (-1+i)/2, (1,1): sp.Rational(3,4)*(1-i),
      (2,1): -sp.Rational(3,4)-sp.Rational(3,4)*i, (1,2): 1+sp.Rational(3,4)*i, (2,2): 1-sp.Rational(3,4)*i }
def tc(side, n, conjg=False):
    if side == 1 and conjg:
        return (-1)**n * t[(2,n)]
    return t[(side,n)]
def rfun(dd, m, conjg=False):
    s = -2*dd if not conjg else 2*dd
    return 2*(dd-2*m) + s*i - 4*(dd-1)*(dd-2*m)/(2*dd-3)
def R2(dd, m, conjg=False):
    s = -1 if not conjg else 1
    return sp.Rational(1,2)*( (3*dd-3)/2 + s*i*rfun(dd,m,conjg)/2 - (dd-2)**2 )
d = 2*p-1
rat = { (p+1,p-1): i*(2*d-1)/(p+1), (p+1,p): (2*d-1)*(2*d+1)/((p+1)*p) }
def e_ratio(key, delta, conjg=False):
    D = key[0]+key[1]
    r0 = rat[key] if not conjg else sp.conjugate(rat[key])
    if delta == 0: return r0
    if delta == 1: return (-1 if not conjg else 1)*i*(D-2)*r0
    if delta == 2: return R2(D, key[0], conjg)*r0
def e_off_slope_ratio(conjg=False):
    return -rat[(p+1,p)]/(2*(d+2)-3) if not conjg else -sp.conjugate(rat[(p+1,p)])/(2*(d+2)-3)
def ehat(delta, s, conjg=False):
    P = p+1
    if (delta,s)==(0,0): return e_ratio((P,p),0,conjg)*tc(1,0,conjg)*P/(4*p-1)
    if (delta,s)==(1,0): return ( e_ratio((P,p),1,conjg)*tc(1,0,conjg)*P/(4*p-2)
                                + e_ratio((P,p),0,conjg)*tc(1,1,conjg)*P/(4*p-1) )
    if (delta,s)==(2,0): return ( e_ratio((P,p),2,conjg)*tc(1,0,conjg)*P/(4*p-3)
                                + e_ratio((P,p),1,conjg)*tc(1,1,conjg)*P/(4*p-2)
                                + e_ratio((P,p),0,conjg)*tc(1,2,conjg)*P/(4*p-1) )
    if (delta,s)==(2,1): return e_off_slope_ratio(conjg)*tc(1,0,conjg)*P/(4*p-3)
def ed(delta, conjg=False):
    if delta in (0,1): return ehat(delta,0,conjg)
    if delta == 2: return ehat(2,0,conjg) - ehat(2,1,conjg)
LHS = R2(d, p) - 1/(2*d-3)
feed_DL1 = (p+1)*( tc(1,2)*rat[(p+1,p-1)]/(2*d-1) + tc(1,1)*e_ratio((p+1,p-1),1)/(2*d-2) + tc(1,0)*e_ratio((p+1,p-1),2)/(2*d-3) )
feed_DL2 = p*( tc(2,2)*ed(0)/(2*d-1) + tc(2,1)*ed(1)/(2*d-2) + tc(2,0)*ed(2)/(2*d-3) )
m = 4*p-6
g_m  = -( ed(0,True)*tc(1,2,True)/(4*p-3) + ed(1,True)*tc(1,1,True)/(4*p-4) + ed(2,True)*tc(1,0,True)/(4*p-5) )
g_m1 = -( ed(0,True)*tc(1,1,True)/(4*p-3) + ed(1,True)*tc(1,0,True)/(4*p-4) )
g_m2 = -( ed(0,True)*tc(1,0,True)/(4*p-3) )
feed_fold = p*( g_m - i*(4*p-7)*g_m1 - sp.binomial(4*p-6,2)*g_m2 )
RHS = feed_DL1 + feed_DL2 + feed_fold
rem = sp.cancel(sp.together(LHS - RHS))
print('E2b-adjacent remainder:', rem)
'''
open('prove_E2b_adjacent.py','w').write(script)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','prove_E2b_adjacent.py'], capture_output=True, text=True, timeout=1500)
print(r.stdout, r.stderr[-2000:] if r.returncode else '')
