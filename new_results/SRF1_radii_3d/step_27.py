
import mpmath as mp, json
mp.mp.dps = 80
h1 = json.load(open('<project>/research/10_frontier/H1_octahedron_end/hi_coeffs_n3_K80.json'))
tau = [mp.mpf(x) for x in h1['eta']]
d   = [mp.mpf(x) for x in h1['d']]
mine = json.load(open('SRF1_n3_K60.json'))
w = [mp.mpf(x) for x in mine['w']]

def recurrence_fit(a, lo, hi):
    """a_{k+2} = P a_{k+1} - Q a_k :  A=sqrt(Q), theta=acos(P/2A)  (conjugate pair model)"""
    S00=S01=S11=T0=T1=mp.mpf(0)
    for k in range(lo, hi+1):
        x0, x1, x2 = a[k], a[k+1], a[k+2]
        S00 += x1*x1; S01 += x1*x0; S11 += x0*x0
        T0 += x1*x2;  T1 += x0*x2
    det = S00*S11 - S01*S01
    P = (T0*S11 - T1*S01)/det
    Q = (S00*T1 - S01*T0)/det
    A = mp.sqrt(Q); th = mp.acos(P/(2*A))
    return A, th, P, Q

print("conjugate-pair recurrence fits: a_{k+2}=P a_{k+1} - Q a_k")
for nm, a, hi in [('tau', tau, 77), ('d', d, 77), ('w', w, 57)]:
    for lo in [30, 45, 60 if hi>=70 else 40]:
        if lo+2 > hi: continue
        A, th, P, Q = recurrence_fit(a, lo, hi)
        print(f"  {nm} k={lo}..{hi}: A={mp.nstr(A,9)}  theta={mp.nstr(th,9)}  t*={mp.nstr(A**-1*mp.e**(1j*th),9) if False else mp.nstr(mp.mpc(mp.cos(th),mp.sin(th))/A,9)}")
tedge = (-mp.log(2)+mp.pi*1j)/2
dedge = tedge/(1-tedge)
print("predicted: tau 1/t*=", mp.nstr(1/abs(tedge),9), "theta=", mp.nstr(mp.arg(tedge),9))
print("predicted: d   1/d*=", mp.nstr(1/abs(dedge),9), "theta=", mp.nstr(mp.arg(dedge),9))
