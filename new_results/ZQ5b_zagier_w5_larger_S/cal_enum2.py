
import sys, time
sys.path.insert(0,'.')
import w5blib as W
F = W.make_field(13, split=[3], inert=[2])
t0=time.time()
cands = W.enum_exceptional_fast(F, B=8)
print(f"D=13 S={{13,2,3}}: {len(cands)} exceptional S-units ({time.time()-t0:.0f}s)")
import collections
print("h dist:", dict(sorted(collections.Counter(c['h'] for c in cands).items())))
for c in cands[:15]:
    print(f"  x={W.fmt_qd(c['x']):>24s} ~ {float(c['x'].emb(30)):+.12f}  vx={c['vx']}  v1m={c['v1m']}  h={c['h']}")
