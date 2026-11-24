import os
os.chdir('<project>/research/08_theory/LEAD2_leading_proof')
txt = open('proof_lead.md').read()
i = txt.find('## 2.'); j = txt.find('## 3.')
print(txt[i:j])
