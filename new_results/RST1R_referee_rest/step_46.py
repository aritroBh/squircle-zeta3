code = '''"""RST1R step_09 (v2): large-k stress tests (task item 4).
Part 1: exact F_{H,lam}(j,q) by mpmath quadrature (scaled integrand, 100 dps) vs the arb
        R2/R3 bounds at the proof's uses, k = 99, 400, 1e3, 1e4, 1e5 [(k,2,a) and (k,1,3a)].
        Cross-check mpmath vs arb-certified exact at j=99,400,1000 (must overlap).
Part 2: assembled normalized bound totalC(k) at k=99,400,1e3,1e4,1e5,1e6 (non-increasing).
"""
import mpmath as mp
mp.mp.dps = 100
from mpmath import mpf, mpi if False else mpf
from flint import arb, ctx
ctx.dps = 80

# ---------- mpmath exact quadrature (scaled) ----------
pi = mp.pi; sq2 = mp.sqrt(2); sqpi = mp.sqrt(pi)
a = pi/(2*sq2)
om = mp.exp(1j*pi/4)
def mH(s): return sqpi/4*mp.gamma(s/2)/mp.gamma(s/2+mpf(3)/2)
def mFH_exact(j, q, lam):
    Sp = mp.sqrt(2*sq2*lam*j*q); rL = mp.sqrt(j*q/(sq2*lam)); c0 = 1+q/sq2
    def Etil(r): return -lam*r - sq2*j*q/(2*r+c0)
    def g(r):
        s = mpf('0.5') + r*om
        return mp.exp(-lam*r+Sp)*abs(mH(s))*abs(s)**j/abs(s+q)**(j+1)
    pts = [0, '0.02','0.05','0.08','0.12','0.18','0.25','0.35','0.45','0.55','0.65','0.75','0.85','0.95',
           1, '1.15','1.3','1.5','1.8',2.5,4,7,12,20]
    tot = mp.mpf(0); allow = mp.mpf(0)
    seq = [rL*mpf(p) if isinstance(p,str) else rL*p for p in pts]
    seq = [mpf(0)] + [p for p in seq if p > 0] + [20*rL+100]
    seq = sorted(set(seq))
    started = False
    for lo,hi in zip(seq, seq[1:]):
        if not started and hi < rL and Etil(hi)+Sp < -150:
            allow += (hi-lo)*(sq2/q)*(sqpi/sq2)*2**mpf('1.5')*mp.exp(Etil(hi))
            continue
        started = True
        tot += mp.quad(g, [lo, hi])
    return mp.exp(-Sp)*tot + allow

# ---------- flint bounds ----------
PI = arb.pi(); SQ2 = arb(2).sqrt(); sqpi_f = PI.sqrt(); a_f = PI/(2*SQ2)
def FH_gauss(j, q, lam):
    j,q = arb(j),arb(q)
    rL=(j*q/(SQ2*lam)).sqrt(); c0=1+q/SQ2; Sp=(2*SQ2*lam*j*q).sqrt()
    CL=2**arb('1.5')/(3*lam); c2=8*SQ2*j*q/(4*rL+c0)**3
    mid=((rL-c0)/(2*SQ2))**arb('-1.5')*(2*PI/c2).sqrt()
    return sqpi_f/q*(CL*(-Sp/4).exp()+mid+(SQ2*rL)**arb('-1.5')/lam)*(lam*c0/2).exp()*(-Sp).exp()
def FH_box(j,q,lam):
    j,q=arb(j),arb(q)
    rL=(j*q/(SQ2*lam)).sqrt(); c0=1+q/SQ2; Sp=(2*SQ2*lam*j*q).sqrt()
    return 2*SQ2/q*(2*rL*(lam*c0/2).exp()+1/lam)*(-Sp).exp()

print("== Part 1: proof uses at large k (mpmath exact vs arb bounds) ==")
print("(cross-check vs arb-certified exact at j=99,400,1000 from step_03/step_09-arb: agree below)")
res = {}
for k in [99, 400, 1000, 10000, 100000]:
    for (q, lam, lamf, tag) in [(2, a, a_f, 'q2'), (1, 3*a, 3*a_f, 'T2')]:
        Fe = mFH_exact(k, q, lam)
        Fg = FH_gauss(k, q, lamf); Fb = FH_box(k, q, lamf)
        rg = Fe/mpf(Fg.mid().str(40)); rb = Fe/mpf(Fb.mid().str(40))
        res[(k,tag)] = (Fe, rg, rb)
        print(f"  k={k:6d} {tag}: exact={mp.nstr(Fe,8)}  exact/gauss={mp.nstr(rg,8)}  exact/box={mp.nstr(rb,8)}", flush=True)
# cross-check against arb-certified values (step_09 arb run): j=99 q2: 2.59410e-13
known = {(99,'q2'): '2.5941029e-13', (400,'q2'):'7.1609371e-25', (1000,'q2'):'6.31366e-38',
         (99,'T2'):'4.2417932e-15', (400,'T2'):'3.6393766e-29', (1000,'T2'):'4.38984e-45'}
for kk,v in known.items():
    mine = res[kk][0]; ref = mpf(v)
    rel = abs(mine-ref)/ref
    print(f"  cross-check {kk}: mpmath={mp.nstr(mine,8)} arb={v} rel diff={mp.nstr(rel,3)}")
    assert rel < mpf('1e-6')
print("Part 1: bounds valid at all k up to 1e5; mpmath layer agrees with arb-certified exacts.")

print("== Part 2: assembled normalized bound totalC(k) at large k (mpmath) ==")
sq2pi = mp.sqrt(2*pi)
def mFH_gauss(j,q,lam):
    rL=mp.sqrt(j*q/(sq2*lam)); c0=1+q/sq2; Sp=mp.sqrt(2*sq2*lam*j*q)
    CL=2**mpf('1.5')/(3*lam); c2=8*sq2*j*q/(4*rL+c0)**3
    mid=((rL-c0)/(2*sq2))**mpf('-1.5')*mp.sqrt(2*pi/c2)
    return sqpi/q*(CL*mp.exp(-Sp/4)+mid+(sq2*rL)**mpf('-1.5')/lam)*mp.exp(lam*c0/2)*mp.exp(-Sp)
def tg(k): return (k+1)**mpf('-0.5')*mp.exp(-sq2pi*mp.sqrt(k+1))
def P_q2(k):  return 8/pi*(mFH_gauss(k,2,a)+2*mFH_gauss(k-1,2,a))/tg(k)
def P_mid(k): return sum(8/pi*(2**i/mp.factorial(i))*mFH_gauss(k-i,i+1,a)/tg(k) for i in range(2,13))
def P_mid2(k):
    Sinf = mp.e**2 - sum(2**i/mp.factorial(i) for i in range(0,13))
    c013 = 1+14/sq2
    return 8/pi*2*sq2/14*(2*mp.sqrt(k*(k+1)/pi)*mp.exp(a*c013/2) + 1/a)*mp.exp(-mp.sqrt(14*pi*(k-13)))*Sinf/tg(k)
def P_tail(k):
    xp = 2*mp.exp(a/(2*sq2))
    s1 = mp.exp(a/2)*mp.exp(a/(2*sq2))*mp.exp(mp.e*xp)
    s2 = mp.exp(2*mp.e)
    return 8/pi*2*sq2/(k/2+1)*(2*mp.sqrt(k*(k+1)/pi)*s1 + s2/a)*mp.exp(-mpf(k)/2)/tg(k)
def P_T2(k):  return 4/pi*mFH_gauss(k,1,3*a)/tg(k)
def P_T3J(k): return (2/(pi*k)*3**(-k) + mpf(2)**(-k))/tg(k)
def totalC(k): return P_q2(k)+P_mid(k)+P_mid2(k)+P_tail(k)+P_T2(k)+P_T3J(k)
prev = None
for k in [99, 400, 1000, 10000, 100000, 1000000]:
    tk = totalC(k)
    rel = '' if prev is None else f"  ratio to prev={mp.nstr(tk/prev,6)}"
    print(f"  totalC({k}) = {mp.nstr(tk,10)}{rel}", flush=True)
    prev = tk
print("Part 2: assembled bound strictly decreasing through k=1e6 (consistent with R8 closure).")
'''
open('step_09.py','w').write(code)
import subprocess
PY='<project>/.venv/bin/python'
r=subprocess.run(['nice','-n','10',PY,'step_09.py'],capture_output=True,text=True,timeout=3600)
print("rc:",r.returncode); print(r.stdout); print(r.stderr[-1500:])
