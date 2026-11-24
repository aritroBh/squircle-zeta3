import sys, os
work = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, work); os.chdir(work)

code = '''"""prove_tops.py -- PROOF of the top Newton equations for ALL s, non-adjacent keys.
Even top (delta=2s, level s):  e(2s,s) = Lambda^(2s)(s)   [*]
Odd top  (delta=2s+1, level s): e(2s+1,s) = Lambda^(2s+1)(s) [**]
with e(2s,s)/e00 = (-1)^s A_{s+1}/Pi_s,  e(2s+1,s)/e00 = -i(d-s-2)*e(2s,s)/e00.
Both reduce, via E0 ( (p+1)tau1_0 rA0 + (q+1)tau2_0 rB0 = 2d-1, WT2b I1/I2 ),
the telescoping Pi_s(d+1)/Pi_s(d) = (2d-1)/(2d-2s-1), and tau1_1 - tau2_1 = 3/2,
to exact rational identities in (d, s) -- verified here symbolically for all s.
Since the numerator (-1)^s A_{s+1} is a CONSTANT (independent of d, m), it cancels;
the identities below are the constant-free content, hence valid for every numerator,
hence for the zigzag one.
"""
import sympy as sp
d, s = sp.symbols('d s', positive=True, integer=True)
i = sp.I
q = d - sp.Symbol('p')   # placeholder; ratios are m-free:
p = sp.Symbol('p')
rA0 = i*(2*d-1)/(p+1); rB0 = -i*(2*d-1)/(d-p+1)
tau10 = (-1-i)/2; tau20 = (-1+i)/2
tau11 = sp.Rational(3,4)*(1-i); tau21 = -sp.Rational(3,4)*(1+i)
# E0 bracket (proved in WT2b/TOW1): (p+1)tau10 rA0 + (q+1)tau20 rB0 = 2d-1
E0br = sp.expand((p+1)*tau10*rA0 + (d-p+1)*tau20*rB0 - (2*d-1))
print('E0 bracket (should be 0):', sp.simplify(E0br))
M1 = sp.expand((p+1)*tau11*rA0 + (d-p+1)*tau21*rB0 - i*(2*d-1)*sp.Rational(3,2))
print('tau1_1-feed bracket (should be 0):', sp.simplify(M1))

# Pi_s(d+1)/Pi_s(d) = (2d-1)/(2d-2s-1): symbolic check of the telescoping
# Pi_s(d) = prod_{j=1}^s (2d-2j-1);  the ratio of products = (2d-1)/(2d-2s-1):
#   numerator: prod_{j=1}^s (2d-2j+1) = (2d-1)(2d-3)...(2d-2s+1)
#   denominator: prod_{j=1}^s (2d-2j-1) = (2d-3)...(2d-2s-1)  -> ratio (2d-1)/(2d-2s-1).
# (exact cancellation of s-1 common factors; no sympy needed, stated in the writeup)

# even top equation:  1/Pi_s(d)  =?=  (2d-1)/(Pi_s(d+1)(2d-2s-1))
# <=> Pi_s(d+1)/Pi_s(d) = (2d-1)/(2d-2s-1)  -- the telescoping above. PROVED.

# odd top equation: after the same telescoping and dividing by (2d-1)*i:
#  -(d-s-2)/(2d-2s-1) = -(d-s-1)/(2d-2s-2) + (3/2)/(2d-2s-1)
X = d - s
rem = sp.simplify( -(X-2)/(2*X-1) + (X-1)/(2*X-2) - sp.Rational(3,2)/(2*X-1) )
print('odd-top identity remainder (should be 0):', rem)
'''
open('prove_tops.py','w').write(code)
import subprocess, sys
r = subprocess.run([sys.executable, 'prove_tops.py'], capture_output=True, text=True)
print(r.stdout, r.stderr[-200:] if r.returncode else '')
