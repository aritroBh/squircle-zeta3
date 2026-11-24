
import mpmath as mp
mp.mp.dps = 100

def load(fn):
    out = {}
    for line in open(fn):
        parts = line.split("\t")
        k = int(parts[0]); mid = parts[1].strip()
        out[k] = mp.mpf(mid.split(" ")[0].strip("[]"))
    return out

d = load("mk_deg9.tsv")
spi = mp.sqrt(mp.pi)

# Natural model:  Y_k := -M_k * N e^theta / sqrt(pi) = cos(theta) [1 + c_r N^{-1/2} + d_r N^{-1} + ...]
#                                                    + sin(theta) [c_i N^{-1/2} + d_i N^{-1} + ...]
# theta = sqrt(pi N), N = k+1.  Predicted: c_r = sqrt(pi)/4*(1+pi/12), c_i = sqrt(pi)/4*(1-pi/12).
cr_pred = spi/4*(1+mp.pi/12); ci_pred = spi/4*(1-mp.pi/12)
print("predicted c_r =", mp.nstr(cr_pred, 20))
print("predicted c_i =", mp.nstr(ci_pred, 20))

Y = {}
for k in range(2, 401):
    N = mp.mpf(k+1); th = mp.sqrt(mp.pi*N)
    Y[k] = (-d[k])*N*mp.exp(th)/spi

def fit(kmin, nterms, kmax=400):
    # unknowns: c_r, c_i, d_r, d_i, e_r, e_i ... pairs at orders N^{-1/2}, N^{-1}, N^{-3/2}
    ks = list(range(kmin, kmax+1))
    m = 2*nterms
    A = mp.matrix(len(ks), m); y = mp.matrix(len(ks),1)
    for i,k in enumerate(ks):
        N = mp.mpf(k+1); th = mp.sqrt(mp.pi*N)
        base = mp.cos(th)
        y[i] = Y[k] - base
        for j in range(nterms):
            powr = N**(-(mp.mpf(j)+mp.mpf('0.5')))
            A[i, 2*j]   = mp.cos(th)*powr
            A[i, 2*j+1] = mp.sin(th)*powr
    G = A.T*A; rhs = A.T*y
    sol = mp.lu_solve(G, rhs)
    rms = mp.mpf(0)
    for i in range(len(ks)):
        pred = sum(A[i,jj]*sol[jj] for jj in range(m))
        rms += (pred-y[i])**2
    return sol, mp.sqrt(rms/len(ks))

for nt in [1,2,3]:
    for kmin in [40, 60, 100, 150, 200]:
        sol, rms = fit(kmin, nt)
        cr, ci = sol[0], sol[1]
        line = f"nterms={nt} kmin={kmin:3d}: c_r={mp.nstr(cr,14)} (err {mp.nstr(cr-cr_pred,3)})  c_i={mp.nstr(ci,14)} (err {mp.nstr(ci-ci_pred,3)})"
        if nt >= 2:
            line += f"  d_r={mp.nstr(sol[2],8)} d_i={mp.nstr(sol[3],8)}"
        line += f"  rms={mp.nstr(rms,2)}"
        print(line)
    print()
