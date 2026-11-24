
import os
BASE = '<project>'
WORK = os.path.join(BASE, 'research/08_theory/WT3_weight3_all_k')
os.chdir(WORK)
with open('progress.log','a') as fh:
    fh.write("step_C: reduction theorem PROVED all k: C1-C4 exact PF identities PASS (affine Gm law => affine kernel law, no +-i poles)\n")
    fh.write("step_04: instrumented reducer k=2..15; W1 (fixed-point deg-2 eqns) PASS all k; W2 (affine Gm law) PASS k=3..15; record_wt3.pkl saved\n")
# Look at WT2b key step files: step_31..34 (tower data + rho discovery), step_39 (recursion)
WT2B = os.path.join(BASE, 'research/08_theory/WT2b_stabilization')
for f in ['step_39.py']:
    s = open(os.path.join(WT2B,f)).read()
    print(f"===== {f} ({len(s)} chars) =====")
    print(s[:7000])
