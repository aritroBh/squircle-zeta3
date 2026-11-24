
import os
BASE = 'research/08_theory'
REV = os.path.join(BASE, 'D13r_lemma_review')
for f in ['step_36.py','step_39.py','step_40.py']:
    print("="*25, f, "="*25)
    print(open(os.path.join(REV,f)).read())
