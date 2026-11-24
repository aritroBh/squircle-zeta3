import os
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
section = open("closed_form_section.md").read()
md = open("result.md").read()
md += """
## TASK 4 RESULT (assembly + local law)

- Full closed form assembled: P^(5)(1)/u = Wt1+Wt2+Wt3+Wt4+Wt5 (Wt0 = 0 exactly).
  At dps 800: closed-form P5/u matches p5_800 to 2.64e-794; P^(5)(1) = u*(that) matches
  p5_800's P5 = 164.1277591133059706510456733124543613051... to 4.88e-794.
- Verifier: verify_oct5d3.py runs in 52.6 s (fresh G-table build, dps 360), ALL CHECKS PASS
  (25 PASS lines). At dps 360 all residuals are 1e-349..1e-363.

### zeta(5), zeta(2)zeta(3) vs prop:order5 (m=4), basis-dependent
  local law:   C_4^(5)   = (33750-23805 sqrt2) cos(pi/8)  = 78.2028433737847736987918987074
               C_4^(2,3) = (-17370+12315 sqrt2) cos(pi/8) = 42.5354327315259866080186809336
  closed form: zeta(5) coeff in P^(5)(1)        = u*(82937921/2944 - 56662239/2944 sqrt2)
                                                    = 1760.7776471224051037379...
               zeta(2)zeta(3) coeff (pi2*z3 folded in) = u*(28661199/2944 + 4602059/2944 sqrt2)
                                                    = 22073.6243104001168144...
  non-local parts (actual - local):  zeta(5):        1682.57480374862033...
                                     zeta(2)zeta(3): 22031.08887766859083...
  These are basis-dependent (the zeta(5) coefficient shifts under the 26-dim atom
  value-kernel, e.g. the 1100-digit representative gives a different value).

### Closure height (51-term relation for R)
  lcm of denominators = 8832 = 2^7 * 3 * 23;  max |numerator| = 281424457 (~2.8e8)
  (LLL integer-relation height 2.2e9 including the target coefficient).
  L(chi8,5) coefficient = -400/3 + 39752/3 sqrt2   (basis-dependent).
  zeta(5) coefficient   = 82937921/2944 - 56662239/2944 sqrt2   (basis-dependent).

""" + section
open("result.md","w").write(md)
print("result.md finalized, %d chars" % len(md))
