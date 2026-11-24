from mpmath import mp, mpf, pi, log, nstr
import json, os
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")
mp.dps = 220
def si2_fast(z):  # direct recurrence, no nsum
    t = z; s = t; k = 0
    while abs(t) > mpf('1e-210'):
        t = -t * z*z * (2*k+1) / ((2*k+3)**2 * (2*k+2))
        s += t; k += 1
    return s
# sanity vs known: si2_series(1)=0.489..., check vs C estimate pieces
print(nstr(si2_fast(mpf(1)), 30))
# compute humps in chunks, checkpoint each chunk
if os.path.exists("humps.json"):
    humps = {int(k): mpf(v) for k,v in json.load(open("humps.json")).items()}
else:
    humps = {}
for m in range(1, 121):
    if m in humps: continue
    humps[m] = si2_fast((m+1)*pi) - si2_fast(m*pi) - (pi/2)*log((m+1)/m)
    if m % 30 == 0:
        json.dump({str(k): mp.nstr(v, 200) for k,v in humps.items()}, open("humps.json","w"))
        print("checkpoint", m)
json.dump({str(k): mp.nstr(v, 200) for k,v in humps.items()}, open("humps.json","w"))
print("done", len(humps), "humps; H_1 =", nstr(humps[1], 20), " H_100 =", nstr(humps[100], 12))
logline("step_17: exact humps H_m, m=1..120, computed via direct recurrence, cached")
