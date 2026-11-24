# step_31: clean CR over quadratic fields (all Fraction), 12-gon & decagon
from fractions import Fraction as Fr
import itertools, pickle

def qsub(x,y,D): return (x[0]-y[0], x[1]-y[1])
def qmul(x,y,D): return (x[0]*y[0]+D*x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def qinv(x,D):
    a,b = x; den = a*a-D*b*b
    assert den != 0
    return (a/den, -b/den)
def qdiv(x,y,D): return qmul(x,qinv(y,D),D)
ONE = (Fr(1),Fr(0)); ZERO = (Fr(0),Fr(0))

def cr4(a,b,c,d,D):
    """cross-ratio (a-c)(b-d)/((a-d)(b-c)); None = infinity."""
    def sub(x,y):
        if x is None: return None
        if y is None: return None
        return qsub(x,y,D)
    n1, n2, d1, d2 = sub(a,c), sub(b,d), sub(a,d), sub(b,c)
    if n1 is None and d1 is None:  # a = inf
        return qdiv(n2, d2, D)
    if n2 is None and d2 is None:  # b = inf
        return qdiv(n1, d1, D)
    if n1 is None and d2 is None:  # c = inf
        return qdiv(n2, d1, D)
    if n2 is None and d1 is None:  # d = inf
        return qdiv(n1, d2, D)
    num = qmul(n1,n2,D); den = qmul(d1,d2,D)
    return qdiv(num,den,D)

def S3orbit(z,D):
    out = {z, qsub(ONE,z,D)}
    if z != ZERO: out.add(qinv(z,D))
    if z != ONE: out.add(qinv(qsub(ONE,z,D),D))
    if z not in (ZERO,ONE):
        out.add(qdiv(z,qsub(z,ONE,D),D))
        out.add(qdiv(qsub(z,ONE,D),z,D))
    return out

def cross_ratios_QF(points, D):
    out = set()
    for quad in itertools.combinations(points, 4):
        val = cr4(*quad, D)
        assert val not in (ZERO, ONE), ("degenerate", quad, val)
        out |= S3orbit(val, D)
    return out

# 12-gon marked points (a + b sqrt3)
p12 = [(Fr(3,2),Fr(1,2)), (Fr(4),Fr(2)), (Fr(-3),Fr(-2)), (Fr(-1,2),Fr(-1,2))]
pts12 = [ZERO, ONE, None] + p12
CR12 = cross_ratios_QF(pts12, 3)
print("12-gon |CR| =", len(CR12))
z8 = [(2,-1),(Fr(-1,2),Fr(1,2)),(0,Fr(1,3)),(7,-4),(1,Fr(-1,2)),(-1,Fr(2,3)),(Fr(1,2),Fr(-1,6)),(Fr(1,3),0)]
z8 = [tuple(Fr(t) for t in z) for z in z8]
names = ["2-s3","(s3-1)/2","s3/3","7-4s3","(2-s3)/2","(2s3-3)/3","(3-s3)/6","1/3"]
for nm,z in zip(names,z8):
    orb = S3orbit(z,3)
    inter = orb & CR12
    dup = []
    for w in CR12:
        w2 = qmul(w,w,3)
        if w2 == z: dup.append(('sqrt',w))
        if (w2[0]==-z[0] and w2[1]==-z[1]): dup.append(('negsqrt',w))
    z2 = qmul(z,z,3)
    if z2 in CR12: dup.append(('square-in-CR',z2))
    print(f"  {nm}: in CR={z in CR12}; orbit cap CR={sorted(inter)}; dup={dup}")

# decagon
p10 = [(Fr(-1,2),Fr(-1,2)), (Fr(3,2),Fr(1,2))]
CR10 = cross_ratios_QF([ZERO,ONE,None]+p10, 5)
print("decagon |CR| =", len(CR10), sorted(CR10))
pickle.dump({"CR12":CR12,"CR10":CR10}, open("CR_regular.pkl","wb"))
