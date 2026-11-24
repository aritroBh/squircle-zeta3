
import mpmath as mp, json
mp.mp.dps = 90
h1 = json.load(open('<project>/research/10_frontier/H1_octahedron_end/hi_coeffs_n3_K80.json'))
tau = [mp.mpf(x) for x in h1['eta']]
d   = [mp.mpf(x) for x in h1['d']]
mine = json.load(open('SRF1_n3_K60.json'))
w = [mp.mpf(x) for x in mine['w']]
A3 = 2/mp.sqrt(mp.pi**2 + mp.log(2)**2)
th_pred = mp.arg((-mp.log(2)+mp.pi*1j)/2)

def est_seq(a):
    b = [a[k]/mp.factorial(k) for k in range(len(a))]
    out = []
    for k in range(10, len(a)-3):
        u = b[k+1]**2 - b[k]*b[k+2]; wv = b[k]*b[k+3] - b[k+1]*b[k+2]
        if u > 0:
            A = u**(1/(2*k)); cth = max(-1,min(1,-wv/(2*A*u)))
            out.append((k, A, mp.acos(cth)))
    return out

def extrap(seq, idx, lo, hi):
    data = [(t[0], t[idx]) for t in seq if lo <= t[0] <= hi]
    S=[[mp.mpf(0)]*3 for _ in range(3)]; T=[mp.mpf(0)]*3
    for k,y in data:
        b=[mp.mpf(1), mp.log(k)/k, 1/k]
        for i in range(3):
            T[i]+=b[i]*y
            for j in range(3): S[i][j]+=b[i]*b[j]
    x = mp.lu_solve(mp.matrix(S), mp.matrix(T))
    return x[0]

print("extrapolation of A_k, theta_k to k=inf (model Y*(1+(b log k + c)/k)):")
for nm, a in [('tau', tau), ('d', d), ('w', w)]:
    seq = est_seq(a)
    hi = max(t[0] for t in seq)
    for lo in [30, 45, 55]:
        A0 = extrap(seq, 1, lo, hi)
        th0 = extrap(seq, 2, lo, hi)
        print(f"  {nm} win {lo}..{hi}: A_inf={mp.nstr(A0,9)}  theta_inf={mp.nstr(th0,9)}")
print(f"  predicted: A={mp.nstr(A3,9)}  theta={mp.nstr(th_pred,9)}")
