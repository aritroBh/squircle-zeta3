
import os
BASE = 'research/08_theory'
REV = os.path.join(BASE, 'D13r_lemma_review')
files = sorted(os.listdir(REV))
print(len(files), "files")
for f in files:
    print(f, os.path.getsize(os.path.join(REV,f)))
