import sys, subprocess
src = r'''
import sys; sys.path.insert(0,'.')
from fractions import Fraction
import json, math
import mpmath as mp
import sympy as sp
import flint
import w5lib as W
import pipe
from w5lib import QD

def unit_only_test(D, Kmax=6, P=830):
    C = W.letters_info(D)
    eps = C['unit']
    rad = C['letters'][0].D
    # candidates x = ± eps^{-k}, k=1..Kmax, keep |x|<1 (all), and record S-unit status of 1-x
    mp.mp.dps = 1200
    args = []
    for k in range(1, Kmax+1):
        for sgn in (1,-1):
            x = eps**(-k)
            if sgn == -1: x = -x
            sx = 1 if sgn==1 else -1
            f = W.factor_sunit(1-x, C)
            args.append((sx, k, x, f is not None))
    names = ["T", "z5"]
    vals = []
    Tv = pipe.load_vals(D, 1150)["T"] if __import__('os').path.exists(pipe.cand_dir(D)+"/vals_1150.pkl") else W.Lchi_hurwitz(D, 1150)/mp.sqrt(rad)
    vals.append(Tv); vals.append(mp.zeta(5))
    for sx, k, x, ok in args:
        names.append(f"Li5({'-' if sx<0 else ''}eps^-{k}){'ok' if ok else 'NOTsun'}")
        vals.append(W.li5_series(x.emb(1200), 1150))
    le = mp.log(eps.emb(1200))
    names += ["le^5", "le^3*pi2", "le*pi2^2"]
    vals += [le**5, le**3*mp.pi**2, le*mp.pi**4]
    print(f"--- D={D} unit-only test (eps={eps}), {len(args)} args, S-unit status:",
          "".join('.' if ok else 'X' for _,_,_,ok in args))
    found = []
    for Pl in (520, 830):
        rows, dt = pipe.lll_rows(vals, Pl)
        rels = []
        for row in rows:
            a = row[:-1]
            mx = max(abs(x) for x in a)
            if mx == 0 or a[0] == 0: continue
            res = pipe.true_residual(a, vals, Pl)
            if abs(res) < mp.mpf('1e-1000'):
                rels.append((mx, a, res))
        rels.sort(key=lambda t: t[0])
        print(f"  P={Pl}: exact relations with a0!=0: {len(rels)}")
        for mx, a, res in rels[:2]:
            print("    a0=", a[0], " max|a|=", mx, " resid=", mp.nstr(res,2))
        found.append(rels)
    if found[0] and found[1]:
        a = found[1][0][1]; a0 = a[0]
        print(f"  YES (D={D}): formula with a0={a0}")
        for i, nm in enumerate(names):
            if a[i]: print(f"    {nm:>22s}  {Fraction(-a[i], a0)}")
    else:
        print(f"  NO exact relation found at either precision (D={D}, searched k=1..{Kmax})")

unit_only_test(8)
unit_only_test(12)
'''
open('steps/step_42_unitonly.py','w').write(src)
r = subprocess.run(['nice','-n','10',sys.executable,'steps/step_42_unitonly.py'],capture_output=True,text=True,timeout=2400)
print(r.stdout); print("STDERR:", r.stderr[-1500:])
with open('progress.log','a') as f: f.write("step_42_unitonly: D=8/D=12 unit-only tests\n")
