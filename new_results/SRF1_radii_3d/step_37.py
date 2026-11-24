
import mpmath as mp, json
mp.mp.dps = 90
h1d = '<project>/research/10_frontier/H1_octahedron_end/'
j2 = json.load(open(h1d+'hi_coeffs_n2_K80.json'))
tau2 = [mp.mpf(x) for x in j2['eta']]
d2   = [mp.mpf(x) for x in j2['d']]
A2 = 2/mp.pi
print("n=2 control: predicted A = 2/pi =", mp.nstr(A2,9), " t* = +-i*pi/2, |t*| =", mp.nstr(mp.pi/2,9))

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
    return mp.lu_solve(mp.matrix(S), mp.matrix(T))[0]

for nm, a in [('tau2', tau2), ('d2', d2)]:
    seq = est_seq(a)
    for lo in [30, 45, 55]:
        A0 = extrap(seq, 1, lo, 76); th0 = extrap(seq, 2, lo, 76)
        print(f"  {nm} win {lo}..76: A_inf={mp.nstr(A0,9)}  theta_inf={mp.nstr(th0,9)}")
print("predicted: A =", mp.nstr(A2,9), " theta = pi/2 =", mp.nstr(mp.pi/2,9))
# root test for radius zero: |a_k|^(1/k)/k -> A/e
print("\nroot test |a_k|^(1/k)/k (should -> A/e =", mp.nstr(A2/mp.e,7), "and A3/e =", mp.nstr(mp.mpf('0.621668')/mp.e,7),"):")
for k in [40, 60, 80]:
    print(f"  n=2 k={k}: {mp.nstr(abs(tau2[k])**(1/k)/k,8)}")
