import pickle, mpmath as mp
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
mp.mp.dps=50
res=pickle.load(open(W+'/gamma_final.pkl','rb'))
print(f"{'a1':>10} {'B/A':>22} {'75/14 dev':>12} {'Z/A':>22} {'-25/3 dev':>12} {'CLASS/A':>12}")
for r in res:
    A=r['z2^2z35']; B=r['z2z37']; Z=r['z2118']; C=r['CLASS']
    r1=B/A; r2=Z/A; r3=C/A
    print(f"{mp.nstr(r['a1'],6):>10} {mp.nstr(r1,18):>22} {mp.nstr(r1-mp.mpf(75)/14,3):>12} {mp.nstr(r2,18):>22} {mp.nstr(r2+mp.mpf(25)/3,3):>12} {mp.nstr(r3,8):>12}")
