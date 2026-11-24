import os, time
lead2 = '<project>/research/08_theory/LEAD2_leading_proof'
os.chdir(lead2)
with open('progress.log','a') as f:
    f.write(f"{time.strftime('%H:%M:%S')} setup: read LEAD1 result.md/genfun.md/verify_lead1.py fully; planned even-d refereeing + odd-d collapse proof via strip GFs\n")
# check the theory folder with kernels exists
THEORY = '<project>/research/08_theory/WTS_weight_k_structure/theory'
print(os.path.exists(THEORY))
print(sorted(os.listdir(THEORY))[:40])
