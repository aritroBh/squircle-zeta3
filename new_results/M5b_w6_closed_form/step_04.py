import os
base='<project>/research/08_theory'
cands=[d for d in os.listdir(base) if 'M5' in d]
print(cands)
print(open(os.path.join(base,'M5a_next_coefficient.brief')).read()[:2000] if os.path.exists(os.path.join(base,'M5a_next_coefficient.brief')) else 'no brief')
