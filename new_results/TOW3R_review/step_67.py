base = '<project>/research/08_theory/'
t = open(base+'WTS_weight_k_structure/theory/result.md').read()
i = t.find('### 2.2')
j = t.find('### 2.4')
print(t[i:j])
