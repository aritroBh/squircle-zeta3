code = r'''
import mpmath as mp
mp.mp.dps = 100

def load(fn):
    out = {}
    for line in open(fn):
        parts = line.split("\t")
        out[int(parts[0])] = mp.mpf(parts[1].strip().split(" ")[0].strip("[]"))
    return out

d = load("mk_deg9.tsv")
spi = mp.sqrt(mp.pi)
br = 1/(2*spi) + spi/4 + mp.pi**mp.mpf('1.5')/48
bi = -1/(2*spi) + spi/4 - mp.pi**mp.mpf('1.5')/48

Y = {}   # Y_k = -M_k N e^th/sqrt(pi) - cos(th)   [= series in 1/sqrt(N) starting at N^{-1/2}]
for k in range(2,401):
    N = mp.mpf(k+1); th = mp.sqrt(mp.pi*N)
    Y[k] = -d[k]*N*mp.exp(th)/spi - mp.cos(th)

def fit(kmin, kmax, nterms, fix_beta):
    # unknowns: pairs (c_r,c_i),(d_r,d_i),(e_r,e_i) at orders N^{-1/2}, N^{-1}, N^{-3/2}
    ks = list(range(kmin, kmax+1))
    m = 2*nterms - (2 if fix_beta else 0)
    A = mp.matrix(len(ks), m); y = mp.matrix(len(ks),1)
    for i,k in enumerate(ks):
        N = mp.mpf(k+1); th = mp.sqrt(mp.pi*N)
        y[i] = Y[k] - (br*mp.cos(th)+bi*mp.sin(th))/mp.sqrt(N) if fix_beta else Y[k]
        col = 0
        for j in range(nterms):
            if fix_beta and j == 0: continue
            pw = N**(-(mp.mpf(j)+mp.mpf('0.5')))
            A[i,col] = mp.cos(th)*pw; A[i,col+1] = mp.sin(th)*pw
            col += 2
    G = A.T*A; rhs = A.T*y
    sol = mp.lu_solve(G, rhs)
    rms = mp.mpf(0)
    for i in range(len(ks)):
        pred = sum(A[i,jj]*sol[jj] for jj in range(m))
        rms += (pred-y[i])**2
    return sol, mp.sqrt(rms/len(ks)), len(ks)-m

print("=== beta FREE (4 or 6 params) vs beta FIXED at analytic (2 or 4 params) ===")
for nt in [2,3]:
    for kmin in [60,100,150,200,250,300]:
        solF, rmsF, dofF = fit(kmin, 400, nt, False)
        solT, rmsT, dofT = fit(kmin, 400, nt, True)
        cr, ci = solF[0], solF[1]
        print(f"nt={nt} kmin={kmin:3d}: free: c_r={mp.nstr(cr,10)} c_i={mp.nstr(ci,10)} rms={mp.nstr(rmsF,3)} | fixed-beta rms={mp.nstr(rmsT,3)}  (dof~{dofT})")
    print()
print("analytic: c_r =", mp.nstr(br,12), " c_i =", mp.nstr(bi,12))
'''
open('step_13.py','w').write(code)
import subprocess, sys
r = subprocess.run([sys.executable,"step_13.py"], capture_output=True, text=True)
print(r.stdout); print(r.stderr[-800:] if r.stderr else "")
