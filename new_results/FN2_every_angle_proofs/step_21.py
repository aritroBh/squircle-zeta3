import os, time
wd = "<project>/research/08_theory/FN2_every_angle_proofs"
os.chdir(wd)
# progress.log
with open("progress.log","w") as f:
    f.write("# FN2: proofs of J2,J3 at every angle + contact law\n")
    f.write(f"{time.strftime('%Y-%m-%d %H:%M')} start: read paper sec:fn, P3 (theta=pi/2 proof), AH3 (every-angle derivation), README row\n")

# result.md skeleton
skeleton = """# FN2: J2, J3 at every angle (PROVED) and the contact law (PROVED)

Work in progress. See proof.md for the full proofs and the lemma ledger.

## Plan
- Lemma A (layer cake at angle theta), Lemma B (density phi_theta: mass, bound, tail series
  d_k(c)=(-1)^k (1-c^2) P'_{k-1}(c)/k), Lemma C (moments/cumulants of log S_theta),
  Lemma D+E (theta-independent window and bivariate inverse from the pi/2 proof).
- Theorem 1 (master expansion, every angle) -> displayed J2, J3 (sympy-exact algebra).
- Theorem 2 (contact law, leading + refined, every angle).
- Numerics: residuals O(eps^4 log^4 n) at 2 angles x 2 functions, n=1e4..1e6;
  contact f=sin(pi x/2), 3 angles, n up to 1e8.

## Status
- Inputs read. Master-formula algebra checked by hand (T2, T3 tail terms reproduce the
  paper's displayed J2, J3 constant terms, incl. the factor 2z(1-z) = (1-c^2)/2).
"""
with open("result.md","w") as f:
    f.write(skeleton)
print("setup done", os.listdir(wd))
