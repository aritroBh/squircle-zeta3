import sys, json, re
sys.dont_write_bytecode = True
sys.path.insert(0, '.')
sys.path.insert(0, '../COR7_classical_corners')
import mpmath as mp
mp.mp.dps = 75
import numeric_cmzv as NC
NC.set_dps(70)
from cor13c_common import (hpl_to_cmzv, li_atom, atom_log2, atom_log3, atom_ipi,
                           atom_zeta2, atom_zeta3, atom_pi2, atom_pi4, MAPS, _fit_terms)
from cor7_lib import H_series
from fractions import Fraction as Fr

def combo(d):
    tot = mp.mpc(0)
    for w, c in d.items():
        tot += mp.mpf(c.numerator) / c.denominator * NC.cmzv_value(w)
    return tot

worst = 0; wword = None; nwords = 0
# B0: the 12 maps send {0,1/3,1,-1,inf} to claimed roots (numeric)
zeta = mp.e ** (mp.pi * mp.j / 3)
def fnum(exps, z):
    # map with f(0)=z^e0, f(1)=z^e1, f(inf)=z^einf, evaluated at z
    e0, eth, e1, em, einf = exps
    a, b, c = zeta**e0, zeta**e1, zeta**einf
    if z == 'inf': return a/b if False else (c*(b-a)) / (b-a)
    return (c*(b-a)*z + a*(c-b)) / ((b-a)*z + (c-b))
maperr = 0
for exps in MAPS:
    e0, eth, e1, em, einf = exps
    for zval, e in [(mp.mpf(1)/3, eth), (-1, em)]:
        got = fnum(exps, zval)
        maperr = max(maperr, abs(got - zeta**e))
print(f"[maps] 12 maps numeric check, worst |f(z)-zeta^e| = {mp.nstr(maperr,3)}")

# B1: HPL words
WORDS = sorted({w for name in ["hpl_fit_k3_(0,0,0).json", "hpl_fit_k3_(1,0,0).json"]
                for e, w, v in _fit_terms(name) if w}, key=lambda x: (len(x), x))
for w in WORDS:
    err = abs(combo(hpl_to_cmzv(w)) - H_series(w, mp.mpf(1)/3))
    nwords += 1
    if err > worst: worst, wword = err, w
print(f"[HPL] {nwords} word conversions at 1/3, worst |resid| = {mp.nstr(worst,3)} at H{list(wword)}")

# B2: Li atoms
LI = [(2, '1/3'), (2, '-1/3'), (3, '1/3'), (3, '-1/2'), (4, '1/3'), (4, '-1/3'),
      (4, '2/3'), (4, '1/2'), (4, '-2')]
worstL = 0; latom = None
for n, z in LI:
    err = abs(combo(li_atom(n, Fr(z))) - mp.polylog(n, mp.mpf(z)))
    if err > worstL: worstL, latom = err, (n, z)
print(f"[Li] {len(LI)} atom conversions, worst |resid| = {mp.nstr(worstL,3)} at Li_{latom[0]}({latom[1]})")

# B3: elementary atoms
checks = [
    ("zeta(2)", combo(atom_zeta2()) - mp.zeta(2)),
    ("zeta(3)", combo(atom_zeta3()) - mp.zeta(3)),
    ("log2", combo(atom_log2()) - mp.log(2)),
    ("log3", combo(atom_log3()) - mp.log(3)),
    ("i*pi", combo(atom_ipi()) - mp.pi * mp.j),
    ("(i*pi)^2 = -pi^2", combo(atom_pi2()) + mp.pi**2),
    ("(i*pi)^4 = pi^4", combo(atom_pi4()) - mp.pi**4),
]
worstA = max(abs(r) for _, r in checks)
for nm, r in checks: print(f"[atom] {nm}: |resid| = {mp.nstr(abs(r),3)}")

# B5: Li4(1/9) distribution input
lhs = mp.polylog(4, mp.mpf(1)/9) - 8*(mp.polylog(4, mp.mpf(1)/3) + mp.polylog(4, -mp.mpf(1)/3))
print(f"[dist] Li4(1/9) - 8(Li4(1/3)+Li4(-1/3)) = {mp.nstr(abs(lhs),3)}")

print("SUMMARY:", mp.nstr(maperr,3), mp.nstr(worst,3), mp.nstr(worstL,3), mp.nstr(worstA,3), mp.nstr(abs(lhs),3))
ok = max(maperr, worst, worstL, worstA, abs(lhs)) < mp.mpf('1e-60')
print("STEP_24", "PASS" if ok else "FAIL")
with open('progress.log','a') as f:
    f.write(f"step_24: 70dps audits maps={mp.nstr(maperr,2)} hpl({nwords})={mp.nstr(worst,2)} li({len(LI)})={mp.nstr(worstL,2)} atoms={mp.nstr(worstA,2)} dist={mp.nstr(abs(lhs),2)} -> {'PASS' if ok else 'FAIL'}\n")
