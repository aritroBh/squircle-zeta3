WD='<project>/research/08_theory/FLD1_field_rule_proof'
s=open(WD+'/result.md').read()
import re
for m in re.finditer(r'^#+ .*', s, re.M):
    print(m.start(), m.group())