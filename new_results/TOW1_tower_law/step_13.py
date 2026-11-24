src = open('../WTS_weight_k_structure/theory/fastred.py').read()
# print the middle part that was truncated
i0 = src.find('def pf_add')
i1 = src.find('def pf_from_sympy')
print(src[i0:i1])
