import mpmath as mp, json
mp.mp.dps = 70
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
data = json.load(open("<project>/research/08_theory/AV5b_layer_theorem/ckpt_newfam.json"))
for cutoff in [1e9, 0.025, 0.013]:
    pts1 = [(mp.mpf(t), mp.mpf(a)-8) for t,a,b in data if mp.mpf(t) <= cutoff]
    An = extract(pts1, 5)
    print(f"cutoff {cutoff}: A4 = {mp.nstr(An[3][0], 30)} resid vs zeta(4)/4: {mp.nstr(An[3][0]-mp.zeta(4)/4, 5)}  A3 resid: {mp.nstr(An[2][0]-2*mp.zeta(3),5)}")
