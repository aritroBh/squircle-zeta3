# step_29: 12-gon & decagon CR sets (exact quadratic arithmetic), alphabet membership
from fractions import Fraction as Fr
import itertools, pickle

def qadd(x,y,D): return (x[0]+y[0], x[1]+y[1])
def qsub(x,y,D): return (x[0]-y[0], x[1]-y[1])
def qmul(x,y,D): return (x[0]*y[0]+D*x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def qinv(x,D):
    a,b = x; den = a*a-D*b*b
    return (a/den, -b/den)
def qdiv(x,y,D): return qmul(x,qinv(y,D),D)
def qneg(x): return (-x[0],-x[1])

def S3orbit(z, D):
    out = set()
    pts = [z, qsub((1,0),z,D)]
    if z != (0,0): pts.append(qinv(z,D))
    if z != (1,0): pts.append(qinv(qsub((1,0),z,D),D))
    if z != (0,0) and z != (1,0):
        pts.append(qdiv(z, qsub(z,(1,0),D), D))
        pts.append(qdiv(qsub(z,(1,0),D), z, D))
    for p in pts: out.add(p)
    return out

def cross_ratios_QF(points, D):
    """points: list of pairs or None (infinity). returns set of pairs."""
    out = set()
    for quad in itertools.combinations(points, 4):
        # one ordering then orbit
        a,b,c,d = quad
        def sub(x,y):
            if x is None: return None  # inf - y = inf
            if y is None: return None
            return qsub(x,y,D)
        def div(x,y):
            if x is None and y is None: return (1,0)
            if x is None: return None
            if y is None: return (0,0)
            if y == (0,0): return None
            return qdiv(x,y,D)
        def mul(x,y):
            if x is None or y is None:
                if x == (0,0) or y == (0,0): return (0,0)
                return None
            return qmul(x,y,D)
        # cr = (a-c)(b-d)/((a-d)(b-c))
        n1, n2 = sub(a,c), sub(b,d)
        d1, d2 = sub(a,d), sub(b,c)
        num = mul(n1,n2); den = mul(d1,d2)
        val = div(num, den)
        if val is None: continue
        out |= S3orbit(val, D)
    return out

# 12-gon
p = [(Fr(3,2),Fr(1,2)), (4,2), (-3,-2), (Fr(-1,2),Fr(-1,2))]
pts = [(0,0),(1,0),None]+p
CR12 = cross_ratios_QF(pts, 3)
print("12-gon |CR| =", len(CR12))
z8 = [(2,-1),(Fr(-1,2),Fr(1,2)),(0,Fr(1,3)),(7,-4),(1,Fr(-1,2)),(-1,Fr(2,3)),(Fr(1,2),Fr(-1,6)),(Fr(1,3),0)]
z8 = [tuple(Fr(t) for t in z) for z in z8]
names = ["2-s3","(s3-1)/2","s3/3","7-4s3","(2-s3)/2","(2s3-3)/3","(3-s3)/6","1/3"]
for nm,z in zip(names,z8):
    print(f"  {nm}: in CR = {z in CR12}")
# decagon
pd = [(Fr(-1,2),Fr(-1,2)), (Fr(3,2),Fr(1,2))]  # -phi, phi^2
CR10 = cross_ratios_QF([(0,0),(1,0),None]+pd, 5)
print("decagon CR =", sorted(CR10))
print("  1/phi in CR:", (Fr(-1,2),Fr(1,2)) in CR10 or None)
phi_pts = {"phi":(Fr(1,2),Fr(1,2)), "-phi":(Fr(-1,2),Fr(-1,2)), "1/phi":(Fr(-1,2),Fr(1,2)), "-1/phi":(Fr(1,2),Fr(-1,2)), "phi^2":(Fr(3,2),Fr(1,2)), "1/phi^2":(Fr(3,2),Fr(-1,2))}
for nm,z in phi_pts.items():
    print(f"  {nm} in CR: {z in CR10}")
import pickle
pickle.dump({"CR12":CR12,"CR10":CR10}, open("CR_regular.pkl","wb"))
