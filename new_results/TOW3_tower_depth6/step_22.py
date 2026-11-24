base = '<project>/research/08_theory'
src = open(f'{base}/WTS_weight_k_structure/theory/result.md').read()
print(len(src))
print(src[:9000])
