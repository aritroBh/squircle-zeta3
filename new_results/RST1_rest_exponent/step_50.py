import os
os.chdir('<project>/research/08_theory/RST1_rest_exponent')
open('progress.log','a').write("proof.md finalized (ledger R0-R11); verify_rst1.py PASS. THEOREM: |M_k-M_main|<=1030 N^-1/2 e^{-sqrt(2 pi N)}, N>=100\n")
# tidy: step_01, step_02 were scaffolding (directory exploration); keep everything else
print(open('progress.log').read())
print("="*60)
for f in sorted(os.listdir('.')):
    print(f, os.path.getsize(f))
