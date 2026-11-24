import sys, os
workdir = '<project>/research/08_theory/TOW2R_review'
import sympy as sp
x = sp.Symbol('x')
# P_r as used in prove_S34.py: read them from the file
src = open(os.path.join(workdir, 'tow2_copy', 'prove_S34.py')).read()
ix = src.find('P = {')
print(src[ix:ix+300])
# define P_r for r<=4 by extracting from file text manually:
P = {}
P[0] = sp.Integer(1)
P[1] = x*(x-1)/2
P[2] = x*(x-1)*(x-2)*(3*x-5)/24
# extract P[3], P[4] from file
import re
m3 = re.search(r'3:\s*(.*?),\s*4:', src[ix:], re.S)
m4 = re.search(r'4:\s*(.*?)\}', src[ix:], re.S)
print('P3 text:', m3.group(1)[:200])
print('P4 text:', m4.group(1)[:200])
