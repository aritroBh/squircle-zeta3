WD='<project>/research/08_theory/FLD1_field_rule_proof'
for f in ['step_01.py','step_02.py','step_03.py']:
    print("="*40, f)
    print(open(WD+'/'+f).read()[:400])
import os
print("transcript size", os.path.getsize(WD+'/transcript.json'))