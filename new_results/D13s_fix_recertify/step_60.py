
import os
ROOT = 'research/08_theory/D13s_fix_recertify'
p = os.path.join(ROOT, 'result.md')
s = open(p).read()
s = s.replace("| H3: identity | step file (mpmath) | worst rel resid $1.77{\\cdot}10^{-31}$ over $t\\in[10^{-10},10^{10}]$ @60 dps |",
              "| H3: identity + boxes | `fixed/verify_mS.py` | resid $1.77{\\cdot}10^{-31}$; upper-box margin $3.93{\\cdot}10^{-21}$ (all $t$); lower-box margin $+1.31{\\cdot}10^{-21}$ on $N\\ge100$ |")
s = s.replace("| H3: upper box, all $t>0$ | step file (mpmath) | $\\min({\\rm box}_{\\rm hi}-{\\rm truth}) = 3.93{\\cdot}10^{-21}$ |\n", "")
s = s.replace("`fixed/verify_TM_final.py`, `fixed/verify_TM_mono.py`,\n`fixed/verify_theorem_truth.py`;",
              "`fixed/verify_TM_final.py`, `fixed/verify_TM_mono.py`,\n`fixed/verify_theorem_truth.py`, `fixed/verify_mS.py`;")
open(p,'w').write(s)
with open(os.path.join(ROOT,'progress.log'),'a') as f:
    f.write("D13s turn14: verify_mS.py added (identity 1.77e-31; upper-box margin 3.9e-21 all t; lower-box margin +1.3e-21 on N>=100). result.md finalized.\n"
            "FINAL: theorem C = 254.2275176 (ranges 194.1267251 / 254.2275176 / 47.0701250). All certificates rerun and reproduced.\n")
print("done")
