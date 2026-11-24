"""r03: the frame-w5 rows that Delta's certificate uses LITERALLY (u = 1, and u = 7 via conjugation): from the
re-run s08 (used_rows_delta_*.pkl in this folder), plus the rows that the u = 3, 5 images of frame row 80 coincide with.
Each is regenerated from its stored provenance (frame_ds.convert = hypothesis (b)) with W5G2R's own domain check (r04)."""
from rc import *
src = open(HERE + '/r04_frames.py').read()
exec(src[src.index('from functools'):src.index("cs = pickle.load")])
LOG = 'r03.log'
V5 = pickle.load(open(T + '/W5C1_certificate_scoping/validated_w5.pkl', 'rb'))
FR = [{tuple(w): c for w, c in r} for r in V5['rows']]
lit = set()
for t in ('pp', 'pm', 'mp', 'mm'):
    for ufam, uidx, *_ in pickle.load(open(HERE + f'/used_rows_delta_{t}.pkl', 'rb')):
        if ufam == 'framew5': lit.add(uidx)
log(f'frame-w5 rows used literally in the sector certificates of Delta: {len(lit)}: {sorted(lit)}', LOG)
key = {tuple(sorted(prim(r).items())): i for i, r in enumerate(FR)}
img80 = {}
for u in (3, 5):
    k = tuple(sorted(prim(actv(FR[80], u)).items())); img80[u] = key.get(k)
log(f'u=3,5 images of frame row 80 equal (up to scalar) validated rows: {img80}', LOG)
OPEN = [4, 5, 7, 9, 10, 12, 14, 15, 16, 17, 52, 53, 54, 56, 58, 59, 61, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 87, 91, 95, 103, 104, 107, 113]
log(f'literal rows among the 37 open rows: {sorted(lit & set(OPEN))}; rows used literally but not open: {sorted(lit - set(OPEN))}', LOG)
check = sorted(lit | {v for v in img80.values() if v is not None})
bad = []; mz = 0
for i in check:
    g, fs = V5['provenance'][i]
    ts, why = domain(g, fs)
    if ts is None: bad.append((i, why)); continue
    if fd.primitive(fd.convert(g, ts, fd.identity(fs))) != V5['rows'][i]: bad.append((i, 'differs'))
    mz = max(mz, max(abs(qf(x)) for f in fs for k, x in f))
log(f'{len(check)} frame-w5 rows regenerated with own domain check: failures {bad}; max |x| {mz:.4f}', LOG)
pickle.dump(dict(lit=sorted(lit), img80=img80, checked=check), open(HERE + '/r03.pkl', 'wb'))
