import os
wts = '<project>/research/08_theory/WTS_weight_k_structure/theory'
with open(os.path.join(wts, 'result.md')) as f:
    content = f.read()
# find Lemma 2 and Lemma 3
import re
ix2 = content.find('Lemma 2')
ix3 = content.find('Lemma 3')
print(content[:3500])
