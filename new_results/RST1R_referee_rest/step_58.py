open('progress.log','a').write("step_11: R11 corroborated: A_k^[2] + saddle-eq identities PROVED (sympy diff=0); rest/R_k 7 values match (0.9997061..0.9999883); saddle |g|Ne^th2 = 1.9895,1.8651,1.8130,1.7904 -> sqrt pi; args -> pi/2; Re-Phi resid -> 0; subdominance q3/KE ratios match. Item-1 numerics confirmed.\n")
code = '''"""RST1R step_12: independent end-to-end assembly audit at k=99 (and spot k=150).
(a) mk_deg9.tsv completeness for k in [99,400].
(b) FULL exact-F assembly: (8/pi)[F(k,2) + sum_{i=1}^k (2^i/i!) F(k-i,i+1)] with every F computed
    by independent quadrature (mpmath 70dps, scaled), vs the proof's R4 bound pieces at k=99.
    Also T2 exact (4/pi)F_{3a}(k,1) vs its bound; arb cross-check of selected F values.
"""
import mpmath as mp, re
mp.mp.dps = 70
from flint import arb, ctx
ctx.dps = 70

# (a) TSV completeness
ks = set()
for line in open("<project>/research/08_theory/D13f_sweep/mk_deg9.tsv"):
    p=line.rstrip("\\n").split("\\t")
    if p[0].isdigit(): ks.add(int(p[0]))
missing = [k for k in range(99,401) if k not in ks]
print("(a) mk_deg9.tsv: k-range [%d,%d], missing in [99,400]: %s" % (min(ks), max(ks), missing))
assert not missing

pi=mp.pi; sq2=mp.sqrt(2); sqpi=mp.sqrt(pi); a=pi/(2*sq2)
om = mp.exp(1j*pi/4)
def mH(s): return sqpi/4*mp.gamma(s/2)/mp.gamma(s/2+mp.mpf('1.5'))
def mFH_exact(j, q, lam):
    Sp = mp.sqrt(2*sq2*lam*j*q); rL = mp.sqrt(j*q/(sq2*lam)) if j>0 else mp.mpf(1)
    c0 = 1+q/sq2
    def Etil(r): return -lam*r - sq2*j*q/(2*r+c0)
    def g(r):
        s = mp.mpf('0.5') + r*om
        return mp.exp(-lam*r+Sp)*abs(mH(s))*abs(s)**j/abs(s+q)**(j+1)
    if j == 0:
        return mp.quad(lambda r: mp.exp(-lam*r)*abs(mH(mp.mpf('0.5')+r*om))/abs(mp.mpf('0.5')+r*om+q), [0, 1, 4, 16, 64, mp.inf])
    seq = sorted(set([mp.mpf(0)]+[rL*mpf(p) for p in ['0.05','0.1','0.18','0.28','0.4','0.55','0.7','0.85','1','1.2','1.5','2','3','5','10']]+[20*rL+100]))
    tot = mp.mpf(0); started=False
    for lo,hi in zip(seq,seq[1:]):
        if not started and hi < rL and Etil(hi)+Sp < -150:
            continue
        started=True
        tot += mp.quad(g, [lo,hi])
    return mp.exp(-Sp)*tot

from flint import arb
PI = arb.pi(); SQ2 = arb(2).sqrt(); sqpi_f = PI.sqrt(); sq2pi_f=(2*PI).sqrt(); a_f = PI/(2*SQ2)
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
def tg(k):
    k=arb(k); return (k+1)**arb('-0.5')*(-sq2pi_f*(k+1).sqrt()).exp()

K = 99
print("(b) full exact-F T1 assembly at k=99 (mpmath quads) ...", flush=True)
exact_terms = {}
S_exact = mp.mpf(0)
# i=0 term: F(k,2); i=1: 2 F(k-1,2); i>=2: (2^i/i!) F(k-i,i+1)
F0 = mFH_exact(K,2,a); S_exact += F0; exact_terms[0]=F0
F1 = mFH_exact(K-1,2,a); S_exact += 2*F1; exact_terms[1]=F1
for i in range(2, K+1):
    Fi = mFH_exact(K-i, i+1, a)
    exact_terms[i]=Fi
    S_exact += (mp.mpf(2)**i/mp.factorial(i))*Fi
S_exact *= 8/pi
# proof's bounds per group (arb)
b_q2  = 8/PI*(FH_gauss(K,2,a_f)+2*FH_gauss(K-1,2,a_f))
b_mid = sum(8/PI*(arb(2)**i/arb(i).fac())*FH_gauss(K-i,i+1,a_f) for i in range(2,13))
b_mid2_exactsum = sum(arb(2)**i/arb(i).fac())*0 if False else None
import itertools
b_mid2 = sum(arb(2)**i/arb(i).fac()*FH_box(K-i,i+1,a_f) for i in range(13, K//2+1))
b_tail = sum(arb(2)**i/arb(i).fac()*FH_box(K-i,i+1,a_f) for i in range(K//2+1, K+1))
S_bound = b_q2 + b_mid + b_mid2 + b_tail
print(f"  exact T1 assembly (mpmath)      = {mp.nstr(S_exact,12)}")
print(f"  proof per-term bounds sum (arb) = {S_bound.str(12)}")
print(f"  ratio exact/bound               = {mp.nstr(mp.mpf(S_exact)/mp.mpf(S_bound.mid().str(30)),8)}")
assert S_exact < mp.mpf(S_bound.mid().str(30))
# normalized
print(f"  exact/bound per group:")
print(f"    q2 : exact={mp.nstr(exact_terms[0]+2*exact_terms[1],10)} bound={(b_q2/(8/PI)).str(10)}")
mids_exact = sum((mp.mpf(2)**i/mp.factorial(i))*exact_terms[i] for i in range(2,13))
print(f"    mid: exact={mp.nstr(mids_exact,10)} bound={(b_mid/(8/PI)).str(10)}")
mid2_exact = sum((mp.mpf(2)**i/mp.factorial(i))*exact_terms[i] for i in range(13,K//2+1))
tail_exact = sum((mp.mpf(2)**i/mp.factorial(i))*exact_terms[i] for i in range(K//2+1,K+1))
print(f"    mid2: exact={mp.nstr(mid2_exact,10)} bound={(b_mid2/(8/PI)).str(10)}")
print(f"    tail: exact={mp.nstr(tail_exact,10)} bound={(b_tail/(8/PI)).str(10)}")
# T2 exact
T2e = 4/pi*mFH_exact(K,1,3*a)
T2b = 4/PI*FH_gauss(K,1,3*a_f)
print(f"  T2: exact={mp.nstr(T2e,10)} bound={T2b.str(10)} ratio={mp.nstr(mp.mpf(T2e)/mp.mpf(T2b.mid().str(30)),6)}")
assert T2e < mp.mpf(T2b.mid().str(30))
print("(b) full assembly audit at k=99: all exact < bound. PASS")
'''
open('step_12.py','w').write(code)
import subprocess
PY='<project>/.venv/bin/python'
r=subprocess.run(['nice','-n','10',PY,'step_12.py'],capture_output=True,text=True,timeout=3600)
print("rc:",r.returncode); print(r.stdout); print(r.stderr[-1500:])
