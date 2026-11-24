wts = '<project>/research/08_theory/WTS_weight_k_structure/theory'
s = open(wts+'/result.md').read()
print(len(s))
import re
i2 = s.find('Lemma 2')
i3 = s.find('Lemma 3')
print(s[max(0,i2-3000):i3+4000])
