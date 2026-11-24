
# step_17.py : consolidation verifications.
import mpmath as mp
mp.mp.dps = 105

def A_series_s(p, s=0, tol=None):
    p = mp.mpf(p); s = mp.mpf(s); tol = tol or mp.mpf('1e-100')
    tot = mp.mpf(0); m = 0
    while True:
        m += 1
        t = mp.gamma((m+s)/p)**2 / (p**2 * mp.gamma(1 + 2*(m+s)/p))
        tot += t
        if m > 30 and abs(t) < tol: break
    return tot

def B_series(p, tol=None):
    p = mp.mpf(p); tol = tol or mp.mpf('1e-100')
    tot = mp.mpf(0); n = -1
    while True:
        n += 1
        t = mp.gamma((2*n+1)/p)**2 / (p**2 * mp.gamma(1 + 2*(2*n+1)/p))
        tot += t
        if n > 30 and abs(t) < tol: break
    return tot

# (a) d9 by residual Richardson
d1 = mp.mpf('1.73300259017088785677509047603076917317936436')
d2 = -mp.zeta(2)/2; d3 = mp.zeta(3)/6
d5 = (mp.zeta(2)*mp.zeta(3)-3*mp.zeta(5))/60
d7 = mp.zeta(7)/14 - mp.zeta(2)*mp.zeta(5)/42 - mp.zeta(2)**2*mp.zeta(3)/140
d9cand = (-428400*mp.zeta(9)+136080*mp.zeta(2)*mp.zeta(7)+79*mp.pi**6*mp.zeta(3)+1134*mp.pi**4*mp.zeta(5)-10080*mp.zeta(3)**3)/1814400
ps = list(range(50, 601, 25))
xs=[]; ys=[]
for p in ps:
    e = mp.mpf(1)/p
    D = mp.zeta(2) - A_series_s(p)
    R = (D - d1*e - d2*e**2 - d3*e**3 - d5*e**5 - d7*e**7)/e**9
    xs.append(e); ys.append(R)
# fit ys = d9 + 0*e + d11*e^2 + ... (even/odd pattern: next terms d11 e^11 -> /e^9 gives e^2)
deg = 4
n=len(xs); M=mp.matrix(n,deg)
for i,x in enumerate(xs):
    for j in range(deg): M[i,j]=x**(2*j)
b=mp.matrix(n,1)
for i,y in enumerate(ys): b[i,0]=y
sol=mp.lu_solve(M.T*M,M.T*b)
print('d9 residual-Richardson:', mp.nstr(sol[0,0],18), ' candidate:', mp.nstr(d9cand,18), ' |diff|', mp.nstr(abs(sol[0,0]-d9cand),3))
print('  (fit also gives d11 ~', mp.nstr(sol[2,0],10), ')')

# (b) B d5 with deg-8 Richardson
psB = list(range(60, 601, 20))
valsB = [(mp.mpf(1)/p, mp.mpf(3)/4*mp.zeta(2) - B_series(p)) for p in psB]
deg=8
n=len(valsB); M=mp.matrix(n,deg)
for i,(x,v) in enumerate(valsB):
    for j in range(deg): M[i,j]=x**(j+1)
b=mp.matrix(n,1)
for i,(x,v) in enumerate(valsB): b[i,0]=v
solB=mp.lu_solve(M.T*M,M.T*b)
print('B deficit coefs d1..d8:', [mp.nstr(solB[j,0],9) for j in range(deg)])
d5B_pred = mp.mpf(7)/60*(3*mp.zeta(5)-mp.zeta(2)*mp.zeta(3))
print('B d5 predicted (7/60)(3z5-z2z3) =', mp.nstr(d5B_pred,12))

# (d) triangle-end c1,c2,c3 independently: fit A(1/(1-tau)) = c0 + c1 t + c2 t^2 + ...
taus = [mp.mpf('0.001')*k for k in range(1,61)]
xs=[]; ys=[]
for t in taus:
    p = 1/(1-t)
    xs.append(t); ys.append(A_series_s(p))
deg=6
n=len(xs); M=mp.matrix(n,deg+1)
for i,x in enumerate(xs):
    for j in range(deg+1): M[i,j]=x**j
b=mp.matrix(n,1)
for i,y in enumerate(ys): b[i,0]=y
solC=mp.lu_solve(M.T*M,M.T*b)
print('triangle fit c0..c6:', [mp.nstr(solC[j,0],12) for j in range(deg+1)])
print('expect: c0=z2/3=0.548311355616, c1=0.623642263729, c2=0.335712456134, c3=0.111148100693')

# (f) A_2(1) formula check
pred = mp.zeta(2)/3 - mp.mpf('0.5') - mp.mpf(1)/24
print('A_2(1) =', mp.nstr(A_series_s(1,s=2),20), ' pred z2/3-13/24 =', mp.nstr(pred,20), ' err', mp.nstr(abs(A_series_s(1,s=2)-pred),3))
