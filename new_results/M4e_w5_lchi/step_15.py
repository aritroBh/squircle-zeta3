import os
p = '<project>/research/08_theory/M4d_w5_nice/w5_search.py'
src = open(p).read()
i = src.find('def evaluate')
j = src.find('def mpl')
k = src.find('def numeric')
l = src.find('def lll')
print("evaluate at", i, "mpl at", j, "numeric at", k, "lll at", l)
print(src[src.find('def eliminate'):6000])
