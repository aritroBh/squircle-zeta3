import subprocess
print(subprocess.run(['cat','<project>/research/08_theory/LEAD2_leading_proof/proof_lead.md'], capture_output=True, text=True).stdout[:9000])
