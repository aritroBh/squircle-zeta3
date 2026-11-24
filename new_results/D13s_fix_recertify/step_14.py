
import os
BASE = 'research/08_theory'
REV = os.path.join(BASE, 'D13r_lemma_review')
for f in ['step_28.py','step_35.py']:
    print("="*25, f, "="*25)
    print(open(os.path.join(REV,f)).read())
