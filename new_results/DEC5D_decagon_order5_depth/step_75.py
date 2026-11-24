
import pickle
import mpmath as mp
mp.mp.dps = 250
parsed = pickle.load(open("parsed_mons_dec.pkl","rb"))
localkeys = set()
for mon, cu in parsed:
    a,b,c,d,e,f = mon
    if c==0 and d==0:
        for n in cu: localkeys.add((n,a,b))
print("local keys:", len(localkeys))
ITAB250 = {}
for k in sorted(localkeys):
    n,a,b = k
    f = lambda x: x**n * mp.log(x)**a * mp.log(1-x)**b
    ITAB250[k] = mp.quad(f, [0,1], method='tanh-sinh')
r5 = mp.sqrt(5)
L2v = mp.log(2); LFv = mp.log((1+r5)/2)
acc = mp.mpf(0)
for mon, cu in parsed:
    a,b,c,d,e,f = mon
    if not (c==0 and d==0) or sum(mon)!=5: continue
    constfac = L2v**e * LFv**f
    for n,(p,q) in cu.items():
        acc += (mp.mpf(str(p)) + mp.mpf(str(q))*r5)*constfac*ITAB250[(n,a,b)]
local250 = 1200*acc
print("local block (250):", mp.nstr(local250, 60))
prev = mp.mpf("-495.63773450055328520236660898040695395892707067482")
print("agrees with 110-digit value:", mp.nstr(abs(local250 - prev), 3))
pickle.dump(str(local250), open("local250.pkl","wb"))
