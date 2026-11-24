import os
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
import mpmath as mp, json
mp.mp.dps = 60
# int_R g3 for F1's g3 = (|s|^3/6)(1-tanh|s|)
g3 = lambda s: (abs(s)**3/6)*(1-mp.tanh(s))
Ig3 = mp.quad(g3, [-mp.inf, 0, mp.inf])
print("int_R g3 =", mp.nstr(Ig3, 30), "  7 zeta(4)/32 =", mp.nstr(7*mp.zeta(4)/32, 30), " resid:", mp.nstr(Ig3-7*mp.zeta(4)/32, 6))
# F1 A4 from the 9-point checkpoint with Neville best-column
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
    An=[]; rem=[(d,v) for d,v in pts]
    for k in range(1,order+1):
        val,err = neville([(d,v/d**k) for d,v in rem]); An.append((val,err)); rem=[(d,v-val*d**k) for d,v in rem]
    return An
data = json.load(open("ckpt_newfam.json"))
pts1 = [(mp.mpf(t), mp.mpf(a)-8) for t,a,b in data]
An = extract(pts1, 5)
print("F1 A4 =", mp.nstr(An[3][0], 30), " err", mp.nstr(An[3][1], 4))
print("zeta(4)/4 =", mp.nstr(mp.zeta(4)/4, 30), " resid:", mp.nstr(An[3][0]-mp.zeta(4)/4, 5))
pts2 = [(mp.mpf(t), mp.mpf(b)-8) for t,a,b in data]
An2 = extract(pts2, 5)
print("F2 A4 =", mp.nstr(An2[3][0], 30), " 7 zeta(4)/8 =", mp.nstr(7*mp.zeta(4)/8, 30), " resid:", mp.nstr(An2[3][0]-7*mp.zeta(4)/8, 5))
