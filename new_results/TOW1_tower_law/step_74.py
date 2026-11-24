verify = r'''"""verify_tower.py — reruns every exact identity used in the TOW1 proof, PASS/FAIL per lemma.
Run: nice -n 10 python verify_tower.py    (target: < 20 min on the loaded machine)
Sections:
  V1  (T1) depth-form recursion vs recorded towers (exact QSI), k <= 13
  V2  initial-data closed forms rho(delta), delta <= 2, vs recorded initial towers, d <= 10
  V3  S(t*), S'(t*) closed forms, d <= 10
  V4  Newton identities E1, E2a, E2b as exact k-free identities, all keys, d <= 11
  V5  T1 for the e-built candidate family at every level (incl. blindness E3 at k=2p), k <= 14
  V6  tower law delta <= 2 vs recorded towers, all keys, k <= 14
  V7  fold identity: child diagonal feed = -p DL2 [Gm(v) - Gm(1/v)], p = 2,3
  V8  affine Gm law with WT3 slopes s30, s21, k = 4..14
  V9  the all-d symbolic proofs (prove_E1_adjacent, prove_E2b_nonadjacent, prove_E2b_adjacent, prove_E3)
"""
import sys, os, subprocess
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from towlib import *

RESULTS = []
def report(name, ok, extra=""):
    RESULTS.append((name, bool(ok), extra))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}  {extra}", flush=True)

# V1: recursion T1 vs data
bad = 0; tested = 0
for k in range(3, 14):
    for d in range(2, k):
        for key in keys_of_degree(d):
            i,j = key
            for delta in range(0, 2*d-3):
                tested += 1
                if u(i,j,k,delta) != Lambda_ij(i,j,k,delta, lambda kk,dp: u(kk[0],kk[1],k,dp)):
                    bad += 1
report("V1 (T1) depth recursion vs towers, k<=13", bad==0, f"{tested} cells, {bad} bad")

# V2: rho closed forms
bad = 0; tested = 0
for d in range(2, 11):
    for key in keys_of_degree(d):
        for delta in (0,1,2):
            if delta > 2*d-4: continue
            tested += 1
            if rho_c(key, delta) != rho(key, delta): bad += 1
report("V2 rho(delta) closed forms, d<=10", bad==0, f"{tested} cells, {bad} bad")

# V3: S(t*), S'(t*)
from genintegrand import S_poly
tstar = QSI(Fr(1,2),0,Fr(1,2),0)
bad = 0
for d in range(1, 11):
    for m in range(0, d+1):
        P = S_poly(d,m)
        val = sum( (QSI(cf)*(tstar**s) for s,cf in enumerate(P)), ZERO )
        Pd = [Fr(s)*cf for s,cf in enumerate(P)][1:]
        vald = sum( (QSI(cf)*(tstar**s) for s,cf in enumerate(Pd)), ZERO )
        if S_tstar(d,m) != val or Sp_tstar(d,m) != vald: bad += 1
report("V3 S(t*), S'(t*) closed forms, d<=10", bad==0, f"{bad} bad")

# V4: Newton identities
bad = check_newton(11)
report("V4 Newton identities E1,E2a,E2b, d<=11", len(bad)==0, f"{len(bad)} bad")

# V5: T1 for e-built family at every level (incl. blindness)
bad = 0
for k in range(3, 15):
    for d in range(2, min(10,k)+1):
        for key in keys_of_degree(d):
            i,j = key
            if k == d: continue
            for delta in range(0, min(2, 2*d-4)+1):
                if U_e(key,k,delta) != Lambda_ij(i,j,k,delta, lambda kk,dp: U_e(kk,k,dp)):
                    bad += 1
report("V5 T1 for e-built family (incl. E3 blindness), k<=14", bad==0, f"{bad} bad")

# V6: tower law vs data
bad = check_tower_law()
report("V6 tower law delta<=2 vs data, k<=14", len(bad)==0, f"{len(bad)} bad")

# V7: fold identity
import sympy as sp
from fastred import pf_to_sympy, pf_mul_DL, pf_fold, pf_scale, DL1, DL2, v as fv
def Gm_actual(key, k):
    H = tower(k)[key]
    hp = pf_zero()
    maxr = max(H)
    hp[1] = [unq4(H.get(r,(0,0,0,0))) for r in range(1, maxr+1)]
    hp[2] = [e.conj_i() for e in hp[1]]
    G = pf_integrate_higher(hp)
    C = pf_eval_at(G, BQ)
    return [G[0], G[1], G[2], -C]
aq, bq = AQ.sympy(), BQ.sympy()
DL2e = 1/(fv-aq) + 1/(fv+bq) - 1/fv
def iszero(e):
    n, dd = sp.fraction(sp.together(e))
    return sp.expand(n) == 0
ok = True
for p in (2,3):
    for k in (2*p, 2*p+1, 2*p+2):
        Gm = Gm_actual((p,p),k)
        direct = pf_scale(pf_mul_DL(Gm, DL2), Fr(-p))
        folded = pf_fold(pf_scale(pf_mul_DL(Gm, DL1), Fr(-p)))
        lhs = pf_to_sympy(pf_add(direct, folded))
        G = pf_to_sympy(Gm)
        rhs = -p*DL2e*(G - G.subs(fv,1/fv))
        ok &= iszero(lhs-rhs)
report("V7 fold identity -p DL2[Gm(v)-Gm(1/v)]", ok)

# V8: affine Gm law
def Gm_from_towerlaw(key, k):
    i, j = key; d = i+j
    ipos = [ZERO]*(2*d-1)
    for delta in range(0, 2*d-3):
        ipos[2*d-2-delta] = U_e(key, k, delta)
    hp = pf_zero()
    hp[1] = ipos[1:]
    hp[2] = [e.conj_i() for e in hp[1]]
    while hp[1] and hp[1][-1].is_zero(): hp[1].pop()
    while hp[2] and hp[2][-1].is_zero(): hp[2].pop()
    G = pf_integrate_higher(hp)
    C = pf_eval_at(G, BQ)
    return [G[0], G[1], G[2], -C]
s30e = -(sp.Rational(1,12))*(1 - 2*sp.sqrt(2)*fv/(fv**2+1))
s30pf = pf_from_sympy(s30e)
s21pf = pf_scale(s30pf, Fr(-3))
ok = True
for k in range(4, 15):
    for key, s in (((3,0),s30pf), ((2,1),s21pf)):
        lhs = pf_sub(Gm_from_towerlaw(key,k), Gm_from_towerlaw(key,3))
        if pf_to_sympy(lhs) != pf_to_sympy(pf_scale(s, Fr(k-3))): ok = False
report("V8 affine Gm law with WT3 slopes, k=4..14", ok)

# V9: the all-d symbolic proofs
for f, name in [("prove_E1_adjacent.py","E1 adjacent"), ("prove_E2b_nonadjacent.py","E2b non-adjacent"),
                ("prove_E2b_adjacent.py","E2b adjacent"), ("prove_E3.py","E3 blindness")]:
    r = subprocess.run([sys.executable, os.path.join(HERE,f)], capture_output=True, text=True, timeout=1500)
    out = r.stdout
    ok = ("remainder: 0" in out) or ("remainder = 0" in out)
    if f == "prove_E3.py":
        ok = out.count("= 0") == 3
    report(f"V9 {name} (symbolic, all d)", ok, out.strip().replace(chr(10),' | '))

print()
nfail = sum(1 for _,ok,_ in RESULTS if not ok)
print(f"verify_tower: {len(RESULTS)-nfail}/{len(RESULTS)} PASS, {nfail} FAIL")
'''
open('verify_tower.py','w').write(verify)
print('verify_tower.py written')
