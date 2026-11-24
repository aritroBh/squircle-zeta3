base = '<project>/research/08_theory'
txt = open(f'{base}/TOW1_tower_law/towlib.py').read()
# print the truncated middle part
i1 = txt.find('def polyval')
i2 = txt.find('def LambdaE_ij')
print(txt[i1:i2])
