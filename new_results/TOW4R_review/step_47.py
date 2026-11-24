
import os, pickle, time
from fractions import Fraction as Fr
base = "<project>"
wts = base + "/research/08_theory/WTS_weight_k_structure/theory"
K = {k: pickle.load(open(f'{wts}/kern_{k}.pkl','rb')) for k in [5,6,7,8,9,10,11,12,13]}
print("kern_7 keys:", list(K[7].keys())[:6])
print("kern_7[(5,0)] fields:", list(K[7][(5,0)].keys()))
print("kern_7[(5,0)]['real'] poles:", list(K[7][(5,0)]['real'].keys()))
print("kern_7[(5,0)]['ipos']:", K[7][(5,0)]['ipos'], " 'ineg':", K[7][(5,0)]['ineg'])
print("sample real pole value:", list(K[7][(5,0)]['real'].items())[0])
