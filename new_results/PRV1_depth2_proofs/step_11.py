
import os
os.chdir('research/08_theory/PRV1_depth2_proofs')
log = open('progress.log','a')
log.write("""
Turn 2-3: Read CX1 verify.py (PASSES), explore.py, BRIEF.md, DEC1 result.md sec.3.
Key structure extracted:
- Octagon: Lam = Li_21(x1, x2/x1), swap = Li_21(x2, x1/x2), x2 = x1/(x1-1) Landen pair.
- G(t) = Li_21(-t, -(1-t)/(2t)) i.e. (x,p)=( -t, (1-t)/2 ), p=xy; H(t) = swap of (x,p) = ((1-t)/2, -t).
- So H(t) = Li_21(p(t), x(t)/p(t)) with the SAME curve parametrization: H = G with x<->p swapped.
- Li_11(x,y) closed form derived: log(1-x)log(1-y) + Li2(-y/(1-y)) - Li2(-y(1-x)/(1-y)).
- Partial derivs: d_x Li_21 = Li_11/x; d_y Li_21 = (Li2(x) - Li2(xy)/y)/(1-y).
DECAGON candidate families (decagon points (-1/phi, 1/phi^2) and (1/phi^2, -1/phi)):
- G_dec(t) = Li_21(-t, -(1-t)/t)  [curve p = 1-t; at t=1/phi: y=-1/phi OK; G_dec(1)=0]
- H_dec(t) = Li_21(1-t, -t/(1-t)) [swap of (x,p); H_dec(0)=0]
Derived candidate closed form for G_dec by hand:
  G_dec(t) = zeta(3) - (1/2) log t * Li2(t^2) + Li3(t^2) - 2 Li3(t) - Li3(1-t)
  (used dilog duplication Li2(t)+Li2(-t)=Li2(t^2)/2; pinned C=zeta(3) at t->1).
  At t=1/phi this gives zeta(3) + l_phi z2/5 - l_phi^3/2 - 2 Li3(1/phi), which equals DEC1's
  LLL form 2z3/5 + 4 l_phi z2/5 - l_phi^3 + 3Li3(-1/phi) + Li3(1/phi) MODULO the sum relation
  Li3(1/phi)+Li3(-1/phi) = (z3 - l_phi z2)/5 + l_phi^3/6 (itself = 3-term Landen + duplication).
Next: numeric check of G_dec closed form; derive H_dec closed form.
""")
log.close()
print("logged")
