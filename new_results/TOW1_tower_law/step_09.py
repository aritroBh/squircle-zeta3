txt = open('../WTS_weight_k_structure/theory/result.md').read()
# print sections 2.4 end and 3
i2 = txt.find('### 2.4')
print(txt[i2:][:7000])
