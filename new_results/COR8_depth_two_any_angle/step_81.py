
import pickle
import mpmath as mp
mp.mp.dps = 50
R = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/fit_resid.pkl','rb'))
B = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/fit_basis.pkl','rb'))
pts = [mp.mpf(x) for x in R['pts']]
rv = [mp.mpf(x) for x in R['rv']]
M = [[mp.mpf(x) for x in col] for col in B['B']]
names = B['names']
n = len(names)
# solve B^T q = rv at the first n points (square), then check at the rest
A = mp.matrix([[M[j][i] for j in range(n)] for i in range(n)])
b = mp.matrix([rv[i] for i in range(n)])
q = mp.lu_solve(A, b)
# residual at all points
maxres = 0
for i in range(len(pts)):
    v = sum(q[j]*M[j][i] for j in range(n)) - rv[i]
    maxres = max(maxres, abs(v))
print("max residual over all 253 points:", mp.nstr(maxres, 3))
# identify rational coefficients
sol = []
for j in range(n):
    qj = q[j]
    if abs(qj) < mp.mpf('1e-30'):
        continue
    r = mp.pslq([qj, mp.mpf(1)], tol=mp.mpf('1e-35'), maxcoeff=10**9)
    sol.append((names[j], mp.nstr(qj, 25), r))
for nm, qs, r in sol:
    print(f"{qs}  {r}  {nm}")
print(len(sol), "nonzero coefficients")
