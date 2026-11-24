import pickle, mpmath as mp
mp.mp.dps=120
res=pickle.load(open('gamma_w12.pkl','rb'))
print("W12 irreducibles along Gamma (raw and normalized):")
print(f"{'a1':>8} {'max|irred|':>14} {'max|c12| scale':>14} {'ratio':>12}")
for r in res:
    a1v=r['a1']; m=r['maxirr']
    # need scale; recompute from irr? scale not saved in gamma_w12.pkl. Use gamma_pt1 for ref
    print(f"{mp.nstr(a1v,8):>8} {mp.nstr(m,10):>14}")
# The scales weren't saved for the batch; recompute ratio using known first point
# Let's just present raw max|irred| and note scale from gamma_pt1
g1=pickle.load(open('gamma_pt1_w12.pkl','rb'))
c12=g1[12]
CLASS=c12[9]+mp.mpf(14)/9*c12[10]
m=max(abs(c12[7]),abs(c12[8]),abs(CLASS),abs(c12[11]))
scale=max(abs(x) for x in c12)
print(f"\nReference point a1=0.01: max|irred|={mp.nstr(m,8)}, scale={mp.nstr(scale,5)}, ratio={mp.nstr(m/scale,6)}")
print("\nproof/extra status:")
import os
print(open('gamma_w12_extra.log').read() if os.path.exists('gamma_w12_extra.log') else "extra not yet")
