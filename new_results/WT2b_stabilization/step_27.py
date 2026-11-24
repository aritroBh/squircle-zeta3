base='<project>/research/08_theory/WT2_weight2_all_k'
src=open(base+'/fastred.py').read()
# print the truncated middle: from pf_add to pf_from_sympy
i0=src.find('def pf_add')
i1=src.find('def pf_from_sympy')
print(src[i0:i1])
