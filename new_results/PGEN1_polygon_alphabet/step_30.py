# step_11: generic octagon, task 1(b): p_j, CR_e, membership of used Li3 args
from fractions import Fraction as Fr
import itertools

def S3orbit_q(z):
    z = Fr(z); pts=[z,1-z]
    if z!=0: pts.append(1/z)
    if z!=1: pts.append(1/(1-z))
    if z!=0 and z!=1: pts += [z/(z-1),(z-1)/z]
    return sorted(set(pts), key=lambda w:(abs(w.numerator),w.denominator,w.numerator))

def cross_ratios_Q(points):
    out=set()
    for a,b,c,d in itertools.permutations(points,4):
        try:
            if a is None: val=(b-d)/(b-c)
            elif b is None: val=(a-c)/(a-d)
            elif c is None: val=(b-d)/(a-d)
            elif d is None: val=(a-c)/(b-c)
            else: val=(a-c)*(b-d)/((a-d)*(b-c))
        except ZeroDivisionError: continue
        out.add(Fr(val))
    return sorted(out, key=lambda w:(abs(w.numerator),w.denominator,w.numerator))

K12 = {0:(Fr(3,2),Fr(-7,12)), 1:(Fr(1,3),Fr(1,9)), 2:(Fr(5,3),Fr(-5,9)), 3:(Fr(1,4),Fr(1,24))}
G8 = [Fr(-9,16),Fr(4,9),Fr(1,5),Fr(1,3),Fr(-5,3),Fr(-1,5),Fr(5,4),Fr(2,5)]
B0 = [Fr(-2,5),Fr(-3,2),Fr(-7,5),Fr(7,12),Fr(1,3),Fr(2,7),Fr(2,3),Fr(-3,5),Fr(5,6),Fr(1,5),
      Fr(1,7),Fr(6,7),Fr(3,7),Fr(3,5),Fr(7,15),Fr(-3,4),Fr(7,8),Fr(10,9),Fr(-7,8),Fr(-2,7)]

def dup_preimages(z, CR):
    """z follows from CR by duplication if z = w^2 or z = -w^2 with w in CR (up to sign):
    Li3(w^2) = 4Li3(w)+4Li3(-w); more generally z = +/- w^2."""
    hits = []
    for w in CR:
        if w*w == z: hits.append(('square', w))
        if -w*w == z: hits.append(('negative square', w))
    return hits

for e in range(4):
    k1,k2 = K12[e]
    p1,p2 = -1/k1, -1/k2
    pts = [Fr(0),Fr(1),None,p1,p2]
    CR = cross_ratios_Q(pts)
    args = B0 if e==0 else G8
    print(f"--- edge {e}: kappa = ({k1},{k2}), p1 = {p1}, p2 = {p2}; |CR_e| = {len(CR)}")
    print("    CR_e =", [str(c) for c in CR])
    for z in args:
        if z in CR:
            print(f"    arg {z}: IN CR_e")
        else:
            hits = dup_preimages(z, CR)
            # also: inversion/reflection are inside CR already; try duplicating CR elements' negatives
            print(f"    arg {z}: not in CR_e; duplication hits: {hits}")
