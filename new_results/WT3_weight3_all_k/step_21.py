
import os
BASE = '<project>'
s = open(os.path.join(BASE, 'paper/squircle_note/scripts/kummer.py')).read()
print(s[:3000])
print("...")
# signatures in w2_integrals
s2 = open(os.path.join(BASE, 'paper/squircle_note/scripts/w2_integrals.py')).read()
import re
for m in re.finditer(r'^def (\w+)\(([^)]*)\)', s2, re.M):
    print("w2:", m.group(1), "(", m.group(2), ")")
print("PE reduce_li2_all sig:")
s3 = open(os.path.join(BASE, 'paper/squircle_note/scripts/polylog_exact.py')).read()
i = s3.find('def reduce_li2_all')
print(s3[i:i+400])
i = s3.find('def li3_closure')
print(s3[i:i+200])
