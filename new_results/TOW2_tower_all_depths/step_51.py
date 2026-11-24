import sys, os
work = '<project>/research/08_theory/TOW2_tower_all_depths'
txt = open(f'{work}/towlib.py').read()
i1 = txt.find('def LambdaE_ij')
i2 = txt.find('def Ediag_conv')
print(txt[i1:i2])
