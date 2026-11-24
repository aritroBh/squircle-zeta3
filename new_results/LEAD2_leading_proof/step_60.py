import os
os.chdir('<project>/research/08_theory/LEAD2_leading_proof')
txt = open('proof_lead.md').read()
i = txt.find('Step 1'); j = txt.find('Step 4')
print(txt[i:j])
