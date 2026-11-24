
import mpmath as mp, json
mp.mp.dps = 90
h1 = json.load(open('<project>/research/10_frontier/H1_octahedron_end/hi_coeffs_n3_K80.json'))
tau = [mp.mpf(x) for x in h1['eta']]
d   = [mp.mpf(x) for x in h1['d']]
mine = json.load(open('SRF1_n3_K60.json'))
w = [mp.mpf(x) for x in mine['w']]
A3 = 2/mp.sqrt(mp.pi**2 + mp.log(2)**2)
Ad_expected = A3  # binomial-transform argument: same rate for d-series

def hankel_est(a, ks):
    """u_k = a_{k+1}^2 - a_k a_{k+2} ~ (k^beta A^k sin theta)^2 : (u_k)^(1/2k) -> A"""
    out = []
    for k in ks:
        u = a[k+1]**2 - a[k]*a[k+2]
        out.append((k, u, u**(1/(2*k)) if u>0 else None))
    return out

print("Hankel u_k sign check + (u_k)^(1/2k):")
for nm, a in [('tau', tau), ('d', d), ('w', w)]:
    est = hankel_est(a, [10,20,30,40,50,60,70,76] if nm!='w' else [10,20,30,40,50,56])
    print(f" {nm}: ", "  ".join(f"k={k}:{'+' if u>0 else '-'}{mp.nstr(v,7) if v else 'neg'}" for k,u,v in est))
print("predicted A3 =", mp.nstr(A3,9))

# refine: fit u_k = A^{2k} sin^2(th) k^{2 beta} Z^2 -> log u_k = 2k log A + 2 beta log k + const
def fit_u(a, lo, hi):
    ks = list(range(lo, hi+1))
    ys = [mp.log(a[k+1]**2 - a[k]*a[k+2]) for k in ks]
    # LS: y = 2 logA * k + 2 beta log k + c
    S=[[mp.mpf(0)]*3 for _ in range(3)]; T=[mp.mpf(0)]*3
    for k,y in zip(ks,ys):
        b=[k, mp.log(k), mp.mpf(1)]
        for i in range(3):
            T[i]+=b[i]*y
            for j in range(3): S[i][j]+=b[i]*b[j]
    x = mp.lu_solve(mp.matrix(S), mp.matrix(T))
    return mp.e**(x[0]/2), x[1]/2
for lo in [30,40,50,60]:
    if lo+2<=76:
        At, bt = fit_u(tau, lo, 76); Ad_, bd_ = fit_u(d, lo, 76)
        print(f"window {lo}..76: tau A={mp.nstr(At,8)}, beta={mp.nstr(bt,6)} | d A={mp.nstr(Ad_,8)}, beta={mp.nstr(bd_,6)}")
Aw, bw = fit_u(w, 35, 55)
print(f"w window 35..55: A={mp.nstr(Aw,8)}, beta={mp.nstr(bw,6)}")
