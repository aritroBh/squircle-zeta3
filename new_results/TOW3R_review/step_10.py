base = '<project>/research/08_theory/'
txt = open(base+'CX8_tow2_review/report.md').read()
print(len(txt))
# print headings
import re
for m in re.finditer(r'^#+ .*', txt, re.M):
    print(m.group(0))
print('-----')
print(txt[:5000])
