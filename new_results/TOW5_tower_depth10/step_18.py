with open('<project>/research/08_theory/WTS_weight_k_structure/theory/result.md') as f:
    txt = f.read()
# find section 2.3
i23 = txt.find('### 2.3')
i24 = txt.find('### 2.4')
print(txt[i23:i24])
