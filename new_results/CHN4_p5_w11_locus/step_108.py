import pickle, mpmath as mp
mp.mp.dps=120
res=pickle.load(open('gamma_final.pkl','rb'))
print("FINAL normalized W12 table along interior Gamma:")
print(f"{'a1':>10} {'z2^2z35':>14} {'z2z37':>14} {'CLASS':>14} {'z2118':>12} {'ratio':>10}")
for r in sorted(res,key=lambda r:float(r['a1'])):
    print(f"{mp.nstr(r['a1'],8):>10} {mp.nstr(r['z2^2z35'],7):>14} {mp.nstr(r['z2z37'],7):>14} {mp.nstr(r['CLASS'],7):>14} {mp.nstr(r['z2118'],7):>12} {mp.nstr(r['ratio'],5):>10}")
ratios=[float(r['ratio']) for r in res]
print(f"\nmin normalized ratio (max|W12 irred|/max|c12|): {min(ratios):.3g} (at a1={mp.nstr(min(res,key=lambda r:float(r['ratio']))['a1'],6)})")
print(f"max control residual (should be ~0): {mp.nstr(max(r['ctrl'] for r in res),3)}")
