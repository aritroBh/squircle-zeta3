import json
base='research/'
hj=json.load(open(base+'10_frontier/H1_octahedron_end/hi_coeffs_n2_K80.json'))
print("n,K:",hj['n'],hj['K'])
print("eta (first 6):",hj['eta'][:6])
print("len eta:",len(hj['eta']))
print("eta[12]:",hj['eta'][12])
# compare with my exact tau^12 numeric value
import sympy as sp
P=json.load(open('tau_reduced_2d.json'))['tau_reduced']
v=sp.N(sp.sympify(P['12']),40)
print("my exact tau^12 numeric:", v)
