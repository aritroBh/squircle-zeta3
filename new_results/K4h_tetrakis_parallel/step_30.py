import json, glob
import mpmath as mp
mp.mp.dps = 90
rows = []
for f in sorted(glob.glob("ckpt_val96/chunk_*.json")):
    r = json.load(open(f))
    rows.append((r['idx'], mp.mpf(r['acc'][3]), mp.mpf(r['acc'][4])))
print("chunk partials k=4 (and k=3):")
for i, a3, a4 in rows[:8]:
    print("  chunk %02d: k3=%s  k4=%s" % (i, mp.nstr(a3,6), mp.nstr(a4,6)))
m3 = max(abs(r[1]) for r in rows); m4 = max(abs(r[2]) for r in rows)
s4 = sum(r[2] for r in rows)
print("max |chunk k3 partial| =", mp.nstr(m3,3), " max |chunk k4 partial| =", mp.nstr(m4,3))
print("sum k4 partials (=/48) =", mp.nstr(s4,6), " 48x =", mp.nstr(48*s4,6))
print("=> k=4 cancellation ratio ~", mp.nstr(m4/abs(s4),3))
