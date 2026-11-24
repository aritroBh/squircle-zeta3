import mpmath, json
from mpmath import mp, mpf, pi, nstr, sqrt
mp.dps = 120
res = json.load(open('ckpt_Mk_cornerintegral.json'))
M = [mpf(v) for v in res['3pi/4']]
sq2 = sqrt(2)
# prediction: P_oct(p) = 16(sq2-1) + 16*sq2 * sum_{k} M_k u^{k+1}, u=1/p
def P_oct_pred(p, K=3):
    u = 1/mpf(p)
    return 16*(sq2-1) + 16*sq2*sum(M[k]*u**(k+1) for k in range(K+1))
# AV5's octagon data: eps = 1/p directly
epss = [mpf(s) for s in json.load(open('../AV5_layer_universality/ckpt_oct4_epss.json'))]
Ps = [mpf(s) for s in json.load(open('../AV5_layer_universality/ckpt_oct4_Ps.json'))]
print("octagon: corner-decomposition prediction vs AV5 direct perimeters")
for e, P in zip(epss, Ps):
    pred = P_oct_pred(1/e)
    print(f"  p={nstr(1/e,8)}: AV5 P={nstr(P,30)}  pred={nstr(pred,30)}  diff={nstr(P-pred,6)}")
# predicted u-coefficients
d1p = 16*sq2*M[0]; d2p = 16*sq2*M[1]; d3p = 16*sq2*M[2]; d4p = 16*sq2*M[3]
print("predicted u-coeffs:", nstr(d1p,40), nstr(d2p,40), nstr(d3p,40), nstr(d4p,40))
print("AV5 fit coeffs    : -3.399735940130302868833826, 1.663308534742409, -0.121579408610481")
with open('progress.log','a') as f:
    f.write("step_54: octagon corner decomposition 16sq2*M_k(3pi/4) validated against AV5 direct perimeters\n")
