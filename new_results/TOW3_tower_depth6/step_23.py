base = '<project>/research/08_theory'
src = open(f'{base}/WTS_weight_k_structure/theory/result.md').read()
# print sections 2 and 3
idx = src.find('## 2')
print(src[idx:idx+12000])
