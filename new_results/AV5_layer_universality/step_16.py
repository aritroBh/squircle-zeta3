import mpmath, pickle
from mpmath import mp, mpf, sqrt, zeta, log, exp, quad, inf, polylog, nstr
mp.dps = 70
g1 = lambda s: log(1+exp(-2*s))
g3 = lambda s: mpf(2)/3*s**3*(mpmath.tanh(s)-1)

def I(a,b,S=mpf(45)):
    return quad(lambda s: s**a * g1(s)**b, [0, 1, 3, 8, 20, S])
def J(a,b,c,S=mpf(45)):  # s^a g1^b g3^c
    return quad(lambda s: s**a * g1(s)**b * g3(s)**c, [0, 1, 3, 8, 20, S])

ints = {}
for (a,b) in [(0,2),(1,2),(0,3),(1,3),(2,2),(0,4),(1,4),(2,3),(0,5),(3,2),(0,6)]:
    ints[(a,b)] = I(a,b)
    print(f"I({a},{b}) = {nstr(ints[(a,b)],50)}")
for (a,b,c) in [(0,1,1),(1,1,1),(2,1,1),(0,2,1),(1,2,1),(0,3,1)]:
    ints[(a,b,c)] = J(a,b,c)
    print(f"J({a},{b},{c}) = {nstr(ints[(a,b,c)],50)}")
pickle.dump({str(k): nstr(v,60) for k,v in ints.items()}, open('ckpt_profile_integrals.pkl','wb'))

# PSLQ I(0,2) against weight-4 basis: z4, z3 ln2, z2 ln2^2, ln2^4, Li4(1/2)
ln2 = log(2)
basis4 = [zeta(4), zeta(3)*ln2, zeta(2)*ln2**2, ln2**4, polylog(4, mpf(1)/2)]
target = ints[(0,2)]
rel = mpmath.pslq(basis4 + [target])
print("PSLQ I(0,2) vs [z4, z3 ln2, z2 ln2^2, ln2^4, Li4(1/2), -target]:")
print(rel)
