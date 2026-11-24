base='<project>/research/08_theory/WT2_weight2_all_k'
src=open(base+'/fastred.py').read()
i0=src.find('# split:')
i1=src.find('# ---------- conversion from sympy')
print(src[i0:i1])
