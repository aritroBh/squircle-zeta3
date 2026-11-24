
import os
BASE = 'research/08_theory'
REV = os.path.join(BASE, 'D13r_lemma_review')
for f in ['step_22.py','step_23.py','step_29.py','step_30.py']:
    print("="*25, f, "="*25)
    print(open(os.path.join(REV,f)).read())
