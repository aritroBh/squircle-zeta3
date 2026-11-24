l2='<project>/research/08_theory/LEAD2_leading_proof'
import os, re
print(os.listdir(l2))
for f in ['proof_lead.md','proof_part2.md']:
    txt=open(l2+'/'+f).read()
    print(f, len(txt))
    for m in re.finditer(r'^#+ .*', txt, re.M):
        print('  ', m.group(0))
