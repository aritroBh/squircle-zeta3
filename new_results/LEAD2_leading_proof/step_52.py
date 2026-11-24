import os
os.chdir('<project>/research/08_theory/LEAD2_leading_proof')
txt = open('proof_lead.md').read()
i = txt.find('**Theorem 4.2'); j = txt.find('**Referee notes')
print(txt[i:j])
i2 = txt.find('### 5.1'); j2 = txt.find('### 5.2')
print(txt[i2:j2])
