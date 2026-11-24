import sys, json
sys.path.insert(0,'.')
from fractions import Fraction as Fr
from math import comb, factorial
from step_02 import E_series, poch_ratio, const, padd, pmul, pscale, ZV

nv, k = 2, 4
Ek = E_series(nv, k)
P4 = Ek[k]

def apply_op(P, dA, dL):
    R = {}
    for (j, iA, iL), v in P.items():
        if iA >= dA and iL >= dL:
            key = (j, iA-dA, iL-dL)
            fac = factorial(iA)//factorial(iA-dA) * factorial(iL)//factorial(iL-dL)
            R[key] = R.get(key, 0) + v*fac
    return R

def op_on_ab(P, ra, sb):
    ops = {}
    for x in range(ra+1):
        for y2 in range(sb+1):
            cA = comb(ra, x) * (-(nv-1))**(ra-x) * (-1)**x
            cB = comb(sb, y2) * nv**(sb-y2)
            key = (ra-x + sb-y2, x + y2)
            ops[key] = ops.get(key, 0) + cA*cB
    tot = {}
    for (dA, dL), cv in ops.items():
        tot = padd(tot, pscale(apply_op(P, dA, dL), cv))
    return tot

# compare square_3 = dA^3 - da^3 - db^3 vs 3 da^2 db + 3 da db^2
s3 = padd(apply_op(P4,3,0), pscale(op_on_ab(P4,3,0), -1))
s3 = padd(s3, pscale(op_on_ab(P4,0,3), -1))
alt = padd(pscale(op_on_ab(P4,2,1), 3), pscale(op_on_ab(P4,1,2), 3))
diff = padd(s3, pscale(alt, -1))
print('square_3 == 3 da^2 db + 3 da db^2 :', len(diff)==0)

# So operator algebra OK. Then the value:
# [z3]g_{4,2} via (C_3/3!) square_3 then Heval:
from step_05 import make_DH, Heval
DH = make_DH(nv)
def Cw(w): return Fr((-1)**(w+1)*factorial(w-1))
Q = pscale(s3, Cw(3)/Fr(6))
val = Heval(nv, Q, DH)
print('via square_3:', val)
Q2 = pscale(alt, Cw(3)/Fr(6))
val2 = Heval(nv, Q2, DH)
print('via 3da2db+3dab2:', val2)
t = json.load(open('mu_exact_N150_K6.json'))
print('expected [z3]mu_4,2:', t['4']['2']['z3'], ' g = -6*mu:', -6*Fr(t['4']['2']['z3']))
# also what does rung1 predict: -[z2]mu_3,2 -> -6*(-[z2]mu_{3,2}) = 6*[z2]mu_3,2
print('rung1 expected g:', 6*Fr(t['3']['2']['z2']))
