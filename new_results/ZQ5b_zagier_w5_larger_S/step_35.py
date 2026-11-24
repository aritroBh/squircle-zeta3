import os
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
s = open('w5blib.py').read()
fast = r'''

# ---------------- fast enumeration (integer A,B,C representation) ----------------
def _qd_to_abc(x):
    """QD x = a + b sqrt(D) = (A + B sqrt(D))/C, integers, C>0."""
    den = math.lcm(x.a.denominator, x.b.denominator)
    A = int(x.a*den); B = int(x.b*den)
    g = math.gcd(math.gcd(abs(A),abs(B)),den)
    return (A//g, B//g, den//g)

def _abc_mul(u, v, D):
    A1,B1,C1 = u; A2,B2,C2 = v
    A = A1*A2 + B1*B2*D
    B = A1*B2 + B1*A2
    C = C1*C2
    g = math.gcd(math.gcd(abs(A),abs(B)),C)
    if g>1: A//=g; B//=g; C//=g
    return (A,B,C)

def _abc_pow_table(g, B, D):
    """g^e for e in [-B,B] as abc tuples."""
    ga,gb,gc = _qd_to_abc(g)
    gi = _qd_to_abc(g.inv())
    tab = {0:(1,0,1)}
    u=(1,0,1)
    for e in range(1,B+1):
        u = _abc_mul(u,(ga,gb,gc),D); tab[e]=u
    u=(1,0,1)
    for e in range(1,B+1):
        u = _abc_mul(u,gi,D); tab[-e]=u
    return tab

def _abc_allowed_norm1m(u, D, primes):
    """Is N(1-x) a rational S-unit (only primes allowed)? x=(A+B r)/C."""
    A,B,C = u
    num = (C-A)*(C-A) - B*B*D
    den = C*C
    if num==0: return False
    num=abs(num)
    for p in primes:
        while num%p==0 and num>0: num//=p
        while den%p==0: den//=p
    return num==1 and den==1

def _abc_to_qd(u, D):
    A,B,C = u
    return QD(Fraction(A,C), Fraction(B,C), D)

def enum_exceptional_fast(F, B=8, time_limit=1500, progress_every=200000):
    import numpy as np
    D=F['D']; gens=F['letters']; d=len(gens)
    primes=F['rat_primes']
    Lf = np.array([float(mp.log(g.emb(60))) if g.emb(60)>0 else float(mp.log(-g.emb(60))) for g in gens])
    tabs = [_abc_pow_table(g,B,D) for g in gens]
    t0=time.time()
    axes=[np.arange(-B,B+1)]*d
    grids=np.meshgrid(*axes,indexing='ij')
    T=np.stack([g.ravel() for g in grids],axis=1).astype(float)
    logabs = T @ Lf
    keep = T[logabs < 1e-9].astype(int)
    keep = keep[np.any(keep!=0,axis=1)]
    out={}
    for ni, exps in enumerate(keep):
        if progress_every and ni%progress_every==0 and ni>0:
            el=time.time()-t0
            print(f"  {ni}/{len(keep)} ({el:.0f}s), found {len(out)}", flush=True)
            if el>time_limit:
                print(f"  enum TIME LIMIT; stopping early"); break
        u=(1,0,1)
        for gi,e in enumerate(exps):
            if e: u=_abc_mul(u,tabs[gi][int(e)],D)
        # |x|<1 exact: x=(A+B r)/C, need |A+B sqrt(D)| < C.  use sqrt via isqrt on scaled ints
        A,Bb,C = u
        # sign of A+B sqrt(D): compare A+B sqrt D with 0 via exact integer test
        # compute floor of B sqrt(D)
        import math as _m
        # |A + B sqrt D| < C  <=>  (A + B sqrt D)^2 < C^2 when ... safer: evaluate with mpmath at 40 digits (borderline impossible for these heights)
        v = A + Bb*float(D**0.5)
        if not (abs(v) < C): continue
        if not _abc_allowed_norm1m(u, D, primes): continue
        x = _abc_to_qd(u, D)
        for sgn in (1,-1):
            xx = x if sgn==1 else -x
            f = factor_sunit(1-xx, F)
            if f is None: continue
            out[xx]=dict(x=xx, vx=[int(e) for e in exps], v1m=f, h=int(max(abs(e) for e in exps)))
    res=sorted(out.values(), key=lambda r:(r['h'], r['vx'], -r['x'].a))
    print(f"  enum done: {len(keep)} tuples checked, {len(res)} exceptional, {time.time()-t0:.0f}s")
    return res
'''
s += fast
open('w5blib.py','w').write(s)
# rerun calibration with fast enum
cal = r'''
import sys, time
sys.path.insert(0,'.')
import w5blib as W
F = W.make_field(13, split=[3], inert=[2])
t0=time.time()
cands = W.enum_exceptional_fast(F, B=8)
print(f"D=13 S={{13,2,3}}: {len(cands)} exceptional S-units ({time.time()-t0:.0f}s)")
import collections
print("h dist:", dict(sorted(collections.Counter(c['h'] for c in cands).items())))
for c in cands[:15]:
    print(f"  x={W.fmt_qd(c['x']):>24s} ~ {float(c['x'].emb(30)):+.12f}  vx={c['vx']}  v1m={c['v1m']}  h={c['h']}")
'''
open('cal_enum2.py','w').write(cal)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','cal_enum2.py'], capture_output=True, text=True, cwd='.', timeout=1800)
print(r.stdout); print("ERR:", r.stderr[-800:])
