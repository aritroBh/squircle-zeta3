
chn1 = '<project>/research/08_theory/CHN1_chain_bodies'
src = open(f'{chn1}/mzv_exact.py').read()
# print the truncated middle: from "def stuffle" up to "class DSTable"
i = src.find('def stuffle')
j = src.find('class DSTable')
print(src[i:j])
