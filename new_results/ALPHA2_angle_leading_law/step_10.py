import subprocess
r = subprocess.run(['cat','<project>/research/08_theory/LEAD2_leading_proof/proof_lead.md'], capture_output=True, text=True).stdout
print("LEN", len(r))
print(r[12800:13500])
print("=====")
print(r[20000:])
