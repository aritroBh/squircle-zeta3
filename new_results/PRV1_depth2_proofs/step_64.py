wd = 'research/08_theory/PRV1_depth2_proofs'
with open(f'{wd}/result.md') as f:
    txt = f.read()

old = """with $G_c(t)\\to0$ as $t\\to1^-$ and $H_c(t)\\to0$ as $t\\to0^+$."""
new = """with $G_c(t)\\to0$ as $t\\to1^-$ and $H_c(t)\\to0$ as $t\\to0^+$.
(For $0<c\\le1$ both defining series converge on all of $0<t<1$:
$|x|,|xy|\\in\\{t,\\,c(1-t)\\}<1$; for $c>1$ the series converge on
$1-1/c<t<1$ and the families are understood by analytic continuation. The proved
cases have $c\\in\\{\\tfrac12,1\\}$.)"""
assert old in txt
txt = txt.replace(old, new)
with open(f'{wd}/result.md','w') as f:
    f.write(txt)
print("patched")

with open(f'{wd}/progress.log','a') as f:
    f.write("""
verify_full.py COMPLETE + PASSING (ALL CHECKS PASSED). Added CHAIN-RULE ASSEMBLY part
(derives dG,dH,dGd,dHd exactly from Part 0 formulas + 20 rational curve identities).
Scope tests: T_c involution exact for all c; c-values 1/2 (octagon), 1 (decagon), 2/3
(hexagon-type pair (-1/2,1/3)); PSLQ (150 dps, 19-element weight-3 basis) on
Li21(-1/2,-2/3): no genuine relation (candidate with 5e7 coeffs frozen at 2.2e-110 when
re-evaluated at 250 dps => lattice artifact) => closure special to c in {1/2, 1}.
result.md WRITTEN: Thm O (octagon, exact match OCT6), Thm D (decagon, exact match DEC1),
Lemma L (line-curve family + scope). Next: final clean run + wrap-up.
""")
print("log updated")
