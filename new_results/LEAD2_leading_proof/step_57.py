import os
os.chdir('<project>/research/08_theory/LEAD2_leading_proof')
txt = open('proof_lead.md').read()
print(txt[:3300])
