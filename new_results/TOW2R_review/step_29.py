import os
wts = '<project>/research/08_theory/WTS_weight_k_structure/theory'
with open(os.path.join(wts, 'result.md')) as f:
    content = f.read()
ix = content.find('### Lemma 2')
if ix < 0: ix = content.find('Lemma 2')
print(content[ix:ix+3000])
