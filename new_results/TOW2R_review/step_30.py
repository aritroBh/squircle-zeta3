import os
wts = '<project>/research/08_theory/WTS_weight_k_structure/theory'
with open(os.path.join(wts, 'result.md')) as f:
    content = f.read()
# search for section headers mentioning Lemma
import re
for mm in re.finditer(r'#+ .*Lemma.*', content):
    print(mm.group(0))
# find "Lemma 3." or "Lemma 3 ("
for pat in ['Lemma 3.', 'Lemma 3 (', '**Lemma 3', 'Lemma 2.', '**Lemma 2', 'Lemma 4']:
    ix = content.find(pat)
    print(pat, '->', ix)
