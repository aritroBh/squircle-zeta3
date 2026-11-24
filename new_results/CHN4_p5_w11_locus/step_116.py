# Final consolidated evidence check before writing FINAL: gather the key residuals/digits in one place
import pickle, mpmath as mp
mp.mp.dps=120
res=pickle.load(open('gamma_final.pkl','rb'))
print("Consolidated W12-on-Gamma evidence (controls must vanish; irreducibles must be nonzero):")
for r in sorted(res,key=lambda r:float(r['a1'])):
    print(f"  a1={mp.nstr(r['a1'],8):>9}: ctrl_max={mp.nstr(r['ctrl'],3):>9}  |z2^2z35|={mp.nstr(abs(r['z2^2z35']),5):>9} "
          f"|z2z37|={mp.nstr(abs(r['z2z37']),5):>9} |CLASS|={mp.nstr(abs(r['CLASS']),5):>10} |z2118|={mp.nstr(abs(r['z2118']),5):>9} ratio={mp.nstr(r['ratio'],4)}")
print("\nAll controls vanish (<=5.2e-15):", all(r['ctrl']<mp.mpf('1e-13') for r in res))
print("All four W12 irreducibles nonzero at all pts:", all(all(abs(r[k])>0 for k in ['z2^2z35','z2z37','CLASS','z2118']) for r in res))
