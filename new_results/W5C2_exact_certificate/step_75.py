import os
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
with open('progress.log','a') as f:
    f.write('2026-10-29 03:05:00 W5C2 start: read W5C1/A3h5 pipeline; unit-tested Echelon variants (sparse heap, dense, batched-rounds, RREF, LeanEchelon) vs flint; LeanEchelon fastest (8000 rows/15s on pp). step_01 v3 launched (all 4 sectors, LeanEchelon, order engine->framew5->frame234, stopcheck D in span).\n')
# clean verify_W5.py dead code in lemma check
v = open('verify_W5.py').read()
old = """    # RHS: for each sector where orbit valid: contribute chi(g_w,s)*chi(g,s) at word g(rep)... unfold:
    acc = {}
    for ssec in SECTORS:
        cm = ctx.colmap[ssec]
        if k not in cm: continue
        v = chi(g_w, ssec)  # fold value at the orbit's sector column
        # unfolded isotypic vector z with z[g*rep] = chi(g,s)*v; then out[w] += (1/4)*z... but z[w]=chi(g_w,s)*v
        # so out[w] contribution = (1/4)*chi(g_w,s)*v; out over orbit: out[g*rep] += (1/4) chi(g,s) v
        for g in range(4):
            gw = ctx.orbits[k][g][0] if False else None
    # simpler: check via phi_tilde on the folded dict
"""
new = ""
assert old in v
open('verify_W5.py','w').write(v.replace(old, new))
print('verify cleaned, progress logged')
import subprocess
print(subprocess.run(['tail','-2','step_01.out'], capture_output=True, text=True).stdout)