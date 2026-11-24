import os
os.chdir('<project>/research/08_theory/LEAD2_leading_proof')
txt = open('proof_lead.md').read()
i = txt.find('### 5.2')
print(txt[i:i+4200])
