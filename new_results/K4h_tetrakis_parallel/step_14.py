import mpmath as mp
from mpmath import mpf
for dps in (50, 80, 270):
    mp.mp.dps = dps
    worst = 0; ok = True
    x = mpf(1)
    for i in range(300):
        x = mp.sin(x + mp.pi/7) + mp.log(x+2) + mp.sqrt(x+1)
        for extra in (3, 6, 10):
            if mpf(mp.nstr(x, dps+extra)) != x:
                if extra == 10: ok = False
                worst = max(worst, extra)
    print("dps", dps, "nstr(dps+10) round-trip:", "EXACT" if ok else "FAIL", "(max extra that failed:", worst, ")")
# how many digits does str() actually print?
mp.mp.dps = 270
x = mp.sqrt(2)+mp.pi
print("str len:", len(str(x).replace('-','').replace('.','').lstrip('0')), "dps:", mp.mp.dps)
