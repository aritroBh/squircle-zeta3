
import pickle
import mpmath as mp
mp.mp.dps = 45
R = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/fit_resid.pkl','rb'))
B = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/fit_basis.pkl','rb'))
rv = [mp.mpf(x) for x in R['rv']]
M = [[mp.mpf(x) for x in col] for col in B['B']]
names = B['names']
npts = len(rv)
A = mp.matrix(npts, len(names))
for j in range(len(names)):
    for i in range(npts):
        A[i,j] = M[j][i]
b = mp.matrix(npts, 1)
for i in range(npts): b[i,0] = rv[i]
# truncated SVD least squares
U, S, V = mp.svd(A)
print("singular values: min/max", mp.nstr(S[-1],3), mp.nstr(S[0],3), " rank(>1e-20):", sum(1 for s in S if s > S[0]*mp.mpf('1e-25')))
r = sum(1 for s in S if s > S[0]*mp.mpf('1e-25'))
# q = V_r S_r^{-1} U_r^T b
Utb = mp.matrix(r, 1)
for k in range(r):
    Utb[k,0] = sum(U[i,k]*b[i,0] for i in range(npts))
q = mp.matrix(len(names), 1)
for j in range(len(names)):
    q[j,0] = sum(V[j,k]*Utb[k,0]/S[k] for k in range(r))
res = [b[i,0] - sum(A[i,j]*q[j,0] for j in range(len(names))) for i in range(npts)]
nrm = mp.sqrt(sum(x**2 for x in res))
print("|residual| with full basis (truncated SVD):", mp.nstr(nrm, 3))
big = [(names[j], q[j,0]) for j in range(len(names)) if abs(q[j,0]) > mp.mpf('1e-3')]
print(len(big), "coeffs > 1e-3")
for nm, v in big[:20]:
    print("  ", mp.nstr(v, 12), nm)
pickle.dump({'q': [str(q[j,0]) for j in range(len(names))], 'names': names}, open('fit_svd.pkl','wb'))
