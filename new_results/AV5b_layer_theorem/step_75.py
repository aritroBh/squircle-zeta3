import os
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
import mpmath as mp, json
mp.mp.dps = 90

def neville(pts):
    n = len(pts); Q = [[mp.mpf(0)]*n for _ in range(n)]
    for i in range(n): Q[i][0] = pts[i][1]
    for j in range(1,n):
        for i in range(n-j):
            xi, xj = pts[i][0], pts[i+j][0]
            Q[i][j] = (xj*Q[i][j-1] - xi*Q[i+1][j-1])/(xj-xi)
    errs = [abs(Q[0][k]-Q[0][k-1]) for k in range(1,n)]
    best = min(range(1,n), key=lambda k: errs[k-1])
    return Q[0][best], errs[best-1]
def extract(pts, order=5):
    An = []; rem = [(d, v) for d, v in pts]
    for k in range(1, order+1):
        val, err = neville([(d, v/d**k) for d, v in rem])
        An.append((val, err)); rem = [(d, v - val*d**k) for d, v in rem]
    return An

data = json.load(open("ckpt_newfam.json"))
pts1 = [(mp.mpf(t), mp.mpf(a)-8) for t,a,b in data]
pts2 = [(mp.mpf(t), mp.mpf(b)-8) for t,a,b in data]
z3 = mp.zeta(3); z4 = mp.zeta(4); pi2 = mp.pi**2
print("F1 (multiplicative softplus):")
An = extract(pts1, 5)
for k,(v,e) in enumerate(An,1): print("  A%d = %s (err %s)" % (k, mp.nstr(v,32), mp.nstr(e,3)))
print("  targets: A2 = pi^2/3 =", mp.nstr(pi2/3,32))
print("           A3 = 2 zeta(3) =", mp.nstr(2*z3,32))
print("           A4 = -11 zeta(4)/12 =", mp.nstr(-11*z4/12,32))
print("F2 (additive vertex-LSE):")
An2 = extract(pts2, 5)
for k,(v,e) in enumerate(An2,1): print("  A%d = %s (err %s)" % (k, mp.nstr(v,32), mp.nstr(e,3)))
