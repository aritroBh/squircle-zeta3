import pickle
IT = pickle.load(open("itab140.pkl","rb"))
print(len(IT), "integrals")
# check some values
import mpmath as mp
mp.mp.dps = 50
bad = 0
for m, v in IT.items():
    try:
        x = mp.mpf(v)
        if not abs(x) < mp.mpf('1e50') or abs(x) == 0 and m[1:]+(0,) != (0,)*6:
            pass
    except Exception:
        bad += 1
        if bad < 5: print("bad:", m, repr(v)[:60])
print("unparseable/inf values:", bad)
# specific: b>0 monomial
for m in [(0,0,1,0,0,0,0),(0,0,5,0,0,0,0),(0,5,0,0,0,0,0),(8,5,5,1,1,1,1)]:
    print(m, IT[m][:70])
# exact: int_0^1 (log(1-t)) dt = -1; int (log t) = -1; int t^8 (log t)^5 (log(1-t))^5... 
