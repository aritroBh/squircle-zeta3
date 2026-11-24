"""Test: sparse Echelon (Markowitz) on pp sector rows in GENERATED order vs length-sorted."""
import sys, os, time, pickle
sys.dont_write_bytecode=True; sys.path.insert(0,'.')
import numpy as np
from w5c2_lib import Echelon, log, P1
from w5c2_fold import Ctx, sectorize_raw, norm, W5C1
N = int(sys.argv[1]) if len(sys.argv)>1 else 12000
order = sys.argv[2] if len(sys.argv)>2 else 'gen'
p = P1; sec = (1,1)
ctx = Ctx()
t0 = time.monotonic()
rows = []
if order == 'gen':
    for fam, fn in [('std','rels_std.pkl'),('lift1','rels_lift1.pkl'),('lift2old','rels_lift2old.pkl'),
                    ('lift3old','rels_lift3old.pkl'),('lift4old','rels_lift4old.pkl')]:
        R = pickle.load(open(os.path.join(W5C1,fn),'rb'))
        for lab, row in R:
            raw = sectorize_raw(row, sec, ctx)
            if raw: rows.append(tuple(sorted(raw.items())))
        log(f'{fam} folded rows={len(rows)} {time.monotonic()-t0:.0f}s')
    # frame234 then framew5
    R = pickle.load(open(os.path.join(W5C1,'rels_frame234.pkl'),'rb'))
    for lab, row in R:
        raw = sectorize_raw(row, sec, ctx)
        if raw: rows.append(tuple(sorted(raw.items())))
    log(f'frame234 folded rows={len(rows)} {time.monotonic()-t0:.0f}s')
else:
    from w5c2_lib import load_sector_rows
    rows = [r for _,_,r in load_sector_rows(sec)]
NC = len(ctx.sec_cols[sec])
rrows = []
for r in rows[:N]:
    rrows.append((np.array([c for c,_ in r], np.int32), np.array([int(a)%p for _,a in r], np.int64)))
log(f'feeding {len(rrows)} rows order={order}')
E = Echelon(NC, p)
E.feed(rrows, label=f'{order}', t0=t0, logevery=2000)
nnzs = sorted(len(x[0]) for x in E.pivots)
log(f'DONE rows={len(rrows)} pivots={E.rank()} nnz median={nnzs[len(nnzs)//2]} p90={nnzs[int(len(nnzs)*0.9)]} max={nnzs[-1]} ops={E.ops/1e9:.2f}e9 time={time.monotonic()-t0:.0f}s')
