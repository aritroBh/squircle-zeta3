
import os, subprocess
BASE = 'research/08_theory'
REV = os.path.join(BASE, 'D13r_lemma_review')
# grep for TL0 and breg/Btrue and h3/mS in review steps
for pat in ['TL0', 'Btrue', 'breg']:
    print("="*20, "grep", pat, "="*20)
    r = subprocess.run(['grep','-l',pat]+[os.path.join(REV,f) for f in os.listdir(REV) if f.endswith('.py')],
                       capture_output=True, text=True)
    print(r.stdout)
