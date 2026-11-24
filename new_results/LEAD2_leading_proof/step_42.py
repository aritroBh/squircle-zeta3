import os
os.chdir('<project>/research/08_theory/LEAD2_leading_proof')
txt = open('proof_lead.md').read()
i = txt.find('## 6.')
print(txt[i:i+3500])
