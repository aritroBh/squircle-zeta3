
import mpmath as mp, json
mp.mp.dps = 90
h1 = json.load(open('<project>/research/10_frontier/H1_octahedron_end/hi_coeffs_n3_K80.json'))
tau = [mp.mpf(x) for x in h1['eta']]
d   = [mp.mpf(x) for x in h1['d']]
mine = json.load(open('SRF1_n3_K60.json'))
w = [mp.mpf(x) for x in mine['w']]
A3 = 2/mp.sqrt(mp.pi**2 + mp.log(2)**2)
tedge = (-mp.log(2)+mp.pi*1j)/2
dedge = tedge/(1-tedge)

def estimators(a):
    b = [a[k]/mp.factorial(k) for k in range(len(a))]
    res = []
    for k in range(2, len(a)-3):
        u = b[k+1]**2 - b[k]*b[k+2]
        wv = b[k]*b[k+3] - b[k+1]*b[k+2]
        if u > 0:
            A = u**(1/(2*k))
            cos_th = -wv/(2*A*u) if u != 0 else mp.nan
            # note cos_th has correction factors O(1/k); clip
            res.append((k, A, cos_th))
    return res

for nm, a, pred_A, pred_edge in [('tau', tau, A3, tedge), ('d', d, A3, tedge), ('w', w, A3, tedge)]:
    est = estimators(a)
    print(f"--- {nm} (pred A={mp.nstr(pred_A,7)}, |t*|={mp.nstr(abs(pred_edge),7)}, theta={mp.nstr(mp.arg(pred_edge),7)}) ---")
    for k, A, cth in est:
        if k in (20,30,40,50,60,70,74):
            th = mp.acos(max(-1,min(1,cth))) if abs(cth)<=1 else mp.nan
            print(f"  k={k}: A={mp.nstr(A,8)}  cos_th={mp.nstr(cth,7)}  theta={mp.nstr(th,7)}")
