from mpmath import mp, mpf, pi, log, nstr
import json, os
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")
mp.dps = 220
def si2_fast(z):
    t = z; s = t; k = 0
    while abs(t) > mpf('1e-210'):
        t = -t * z*z * (2*k+1)**2 / ((2*k+3)**3 * (2*k+2))
        s += t; k += 1
    return s
print("Si2(1) =", nstr(si2_fast(mpf(1)), 30), " (expect 0.9818...)")
humps = {}
for m in range(1, 121):
    humps[m] = si2_fast((m+1)*pi) - si2_fast(m*pi) - (pi/2)*log((m+1)/m)
json.dump({str(k): mp.nstr(v, 200) for k,v in humps.items()}, open("humps.json","w"))
print("H_1 =", nstr(humps[1], 20))
print("H_100 =", nstr(humps[100], 12), " expect ~ -6/(100pi)^3 =", nstr(-6/(100*pi)**3, 6))
logline("step_18: recurrence fixed; humps recomputed m=1..120")
