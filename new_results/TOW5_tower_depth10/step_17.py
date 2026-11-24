with open('<project>/research/08_theory/WTS_weight_k_structure/theory/result.md') as f:
    txt = f.read()
import re
# Find Lemma 2 and Lemma 3 statements
for m in re.finditer(r'Lemma [23][^\n]*', txt):
    start = max(0, m.start()-200)
    print(txt[start:m.end()+800])
    print('-'*100)
