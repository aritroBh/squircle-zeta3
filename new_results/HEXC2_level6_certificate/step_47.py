import sys, time, json
sys.dont_write_bytecode = True
sys.path.insert(0, '<project>/research/08_theory/HEXC1_hexagon_certificate')
sys.path.insert(0, '.')
import cmzv6
from hexc2_common import target2, target3
import exact_linalg as XL

t0 = time.monotonic()
lower = {}
cert = {}
for n in (1,2,3):
    rows, _ = cmzv6.make_relations(n, lower)
    lower[n] = rows
    cols = cmzv6.words(n)
    if n < 2: continue
    T = {2: target2, 3: target3}[n]()
    prow, r = XL.select_independent_rows(rows, cols)
    sel = [rows[i] for i in prow]
    pcols = XL.select_pivot_cols(sel, cols)
    y = XL.exact_solve(sel, pcols, cols, T)
    ok, bad = XL.verify_combination(sel, y, cols, T)
    print(f'weight {n}: words={len(cols)} rank={r} pivot_cols={len(pcols)} target_support={len(T)} exact_verify={ok} ({time.monotonic()-t0:.1f}s)')
    assert ok, bad
    cert[n] = dict(cols=[list(w) for w in cols], prow=prow,
                   y=[(str(c.numerator), str(c.denominator)) for c in y])
    # pi-related info for T2: coefficient of (-1,0) [zeta2] and (1,1) [from pi^2]
json.dump({str(k): v for k,v in cert.items()}, open('cert_w23.json','w'))
with open('progress.log','a') as f:
    f.write(f'step_45: EXACT certificates weight 2 and 3 verified (targets in span of proved rows)\n')
print('saved cert_w23.json')
