import os
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
import mpmath as mp, json
mp.mp.dps = 80
exec(open("av5b_fit_true.py").read().split("z3 =")[0].replace("ckpt_true_m3.json","ckpt_dPloc_m3.json"))
z3 = mp.zeta(3); pi2 = mp.pi**2; ln2 = mp.log(2)
pts = load("ckpt_dPloc_m3.json")
An = extract(pts, 4)
print("m=3 local model:")
for k,(v,e) in enumerate(An,1): print("  A%d = %s (err %s)" % (k, mp.nstr(v,28), mp.nstr(e,3)))
A3 = An[2][0]
basis = [z3, z3*ln2, pi2*ln2, ln2**3, pi2, ln2**2, ln2, mp.mpf(1)]
M = [A3/mp.sqrt(3)] + basis
rel = mp.pslq(M, tol=mp.mpf('1e-20'), maxcoeff=10**10, maxsteps=10000)
print("PSLQ A3/sqrt3:", rel)
if rel: print(" resid:", mp.nstr(mp.fsum(mp.mpf(r)*v for r,v in zip(rel,M)),6))
print("C3*z3/6 target:", mp.nstr(-9*mp.sqrt(3)/8*z3/6, 28))
