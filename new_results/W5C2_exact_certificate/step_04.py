"""W5C2 step_04: exact verification and certificate build, one sector at a time.
(a) lemma: u = sum_s (1/4) sum_g chi_s(g) g (Phi_tilde_s(u)) checked on ALL 46656 unit vectors;
(b) per-comp Fraction check sum_j y_j R_j = Phi_s(D) on all sector columns;
(c) preimage check: each used sector row = (1/lam)-fold of its unfolded source row (recomputed);
(d) writes certificate_{sec}.json. Usage: step_04.py pp|pm|mp|mm"""
import sys, os, time, pickle, json
sys.dont_write_bytecode=True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
from w5c2_fold import Ctx, sectorize_raw, norm, fold_D, chi, gapply, W5C1, HERE
from w5c2_lib import log

SECNAME = {'pp':(1,1),'pm':(1,-1),'mp':(-1,1),'mm':(-1,-1)}
secname = sys.argv[1]; sec = SECNAME[secname]
t0 = time.monotonic()
ctx = Ctx()
# ---------- (a) lemma on unit vectors ----------
IDS = ctx.IDS; N = ctx.N
SECTORS = [(1,1),(1,-1),(-1,1),(-1,-1)]
def phi_tilde(u_fold_map):
    """u_fold_map: sec -> dict(k->val). Returns unfolded dict word_id->Fraction via
    w -> sum_s (1/4) sum_g chi_s(g) * val placed at g(rep)... build directly:
    z_s[w] = chi_s(g_w) * fold_s[orbit(w)] for valid orbits; out[w] = sum_s (1/4)*... """
    out = {}
    wm_o, wm_g = ctx.wm_o, ctx.wm_g
    colmap = ctx.colmap
    for w in range(N):
        k = wm_o[w]; g_w = wm_g[w]
        acc = F(0)
        for ssec in SECTORS:
            cm = colmap[ssec]
            if k not in cm: continue
            v = u_fold_map[ssec].get(cm[k])
            if v is None: continue
            acc += chi(g_w, ssec) * v
        if acc: out[w] = acc / 4
    return out
def fold_all(u):
    """u: dict word_id->number -> sec -> dict(col->number) (raw fold)."""
    res = {}
    wm_o, wm_g = ctx.wm_o, ctx.wm_g
    for ssec in SECTORS:
        cm = ctx.colmap[ssec]
        d = {}
        for c, a in u.items():
            k = wm_o[c]
            if k not in cm: continue
            kk = cm[k]
            d[kk] = d.get(kk, 0) + chi(wm_g[c], ssec)*a
        res[ssec] = {k2: v for k2, v in d.items() if v}
    return res
# sparse lemma check: for unit vector e_w, unfold(fold(e_w)) == e_w.
# orbit reps and inverse colmap
reps = [min(x[0] for x in ow) for ow in ctx.orbits]
invc = {ssec: list(ctx.sec_cols[ssec]) for ssec in SECTORS}  # inv[sec][sector_col] = orbit id
bad = 0
for w in range(N):
    k = ctx.wm_o[w]; g_w = ctx.wm_g[w]
    back = {}
    for ssec in SECTORS:
        cm = ctx.colmap[ssec]
        if k not in cm: continue
        v = chi(g_w, ssec)   # fold value at sector col cm[k]
        rep = reps[k]
        for x in (rep, ctx.ti[rep], ctx.si[rep], ctx.ti[ctx.si[rep]]):
            back[x] = back.get(x, 0) + chi(ctx.wm_g[x], ssec)*v
    back = {x: F(a, 4) for x, a in back.items() if a}
    if back != {w: F(1)}: bad += 1
    if bad > 3: break
log(f'step_04 {secname}: (a) fold/unfold lemma checked on all {N} unit vectors, bad={bad} ({time.monotonic()-t0:.0f}s)')
assert bad == 0
# ---------- load sector data ----------
rows = pickle.load(open(os.path.join(HERE,f'work/feed_{secname}.pkl'),'rb'))
d2 = pickle.load(open(os.path.join(HERE,f'work/step02_{secname}.pkl'),'rb'))
S, C = d2['S'], d2['C']
ye = pickle.load(open(os.path.join(HERE,f'work/yexact_{secname}.pkl'),'rb'))
fmap = pickle.load(open(os.path.join(HERE,f'work/foldmap_{secname}.pkl'),'rb'))
# unfolded source rows
ENG = {'std':'rels_std.pkl','lift1':'rels_lift1.pkl','lift2old':'rels_lift2old.pkl',
       'lift3old':'rels_lift3old.pkl','lift4old':'rels_lift4old.pkl','frame234':'rels_frame234.pkl'}
urows = {}
def get_urow(ufam, uidx):
    if ufam not in urows:
        urows[ufam] = pickle.load(open(os.path.join(W5C1, ENG[ufam]),'rb'))
    return urows[ufam]
# D unfolded
Dv = {}
for comp in ('rat','sq'):
    D = json.load(open(os.path.join(W5C1,f'D_{comp}.json')))
    dv = {}
    for ws, val in D.items():
        w = tuple(int(x) for x in ws.split(','))
        dv[ctx.IDS[w]] = F(val)
    Dv[comp] = dv
    globals()[f'fold_{comp}'] = fold_D(dv, sec, ctx)
# ---------- (b)(c) per-comp exact check + preimage check ----------
used = {'rat': {}, 'sq': {}}   # j -> y
famcount = {}
preim_cache = {}
mis = 0
for j, si in enumerate(S):
    yr, ys = ye['y_rat'][j], ye['y_sq'][j]
    if not (yr or ys): continue
    fam, famidx, ttuple = rows[si]
    ufam, uidx, ulab, lam = fmap[fam][famidx]
    # (c) recompute fold of the preimage row and compare
    key = (fam, famidx)
    chk = preim_cache.get(key)
    if chk is None:
        if ufam == 'framew5':
            vw = pickle.load(open(os.path.join(W5C1,'validated_w5.pkl'),'rb'))['rows']
            row_list = [(ctx.IDS[w], c) for w, c in vw[uidx]]
            raw = sectorize_raw(row_list, sec, ctx)
            chk = (tuple(sorted(raw.items())), 1)
        else:
            lab, urow = get_urow(ufam, uidx)[uidx]
            raw = sectorize_raw(urow, sec, ctx)
            nt, lam2 = norm(raw)
            chk = (nt, lam2)
            assert lam2 == lam, (key, lam2, lam)
        preim_cache[key] = chk
    stored_tuple = tuple(ttuple)
    if ufam == 'framew5':
        assert chk[0] == stored_tuple, ('framew5 mismatch', key)
    else:
        # stored_tuple = raw/lam i.e. raw = lam*stored_tuple
        assert chk[0] == stored_tuple, ('engine/frame mismatch', key)
    used['rat'][j] = yr; used['sq'][j] = ys
    famcount[ufam] = famcount.get(ufam, 0) + 1
log(f'step_04 {secname}: (c) preimage checks done, used_rows={len(used["rat"])} families={famcount} ({time.monotonic()-t0:.0f}s)')
for comp in ('rat','sq'):
    acc = {}
    tgt = globals()[f'fold_{comp}']
    yv = ye['y_'+comp]
    for j in used[comp]:
        a = yv[j]
        if not a: continue
        for col, val in rows[S[j]][2]:
            acc[col] = acc.get(col, F(0)) + a*val
    resid = 0; maxr = 0
    for col in set(acc) | set(tgt):
        dd = acc.get(col, F(0)) - tgt.get(col, F(0))
        if dd: resid += 1
    log(f'step_04 {secname}: (b) EXACT sector check {comp}: cols={len(tgt)} residual_nonzero={resid} ({time.monotonic()-t0:.0f}s)')
    assert resid == 0
# ---------- write certificate ----------
cert_rows = []
for j in sorted(used['rat'].keys() | used['sq'].keys()):
    fam, famidx, ttuple = rows[S[j]]
    ufam, uidx, ulab, lam = fmap[fam][famidx]
    cert_rows.append({'j': j, 'S_index': int(S[j]), 'fam': fam, 'famidx': int(famidx),
                      'tuple': [[int(c), int(v)] for c, v in ttuple],
                      'y_rat': str(used['rat'].get(j, F(0))), 'y_sq': str(used['sq'].get(j, F(0))),
                      'src': {'ufam': ufam, 'uidx': int(uidx), 'label': str(ulab), 'lam': str(lam)}})
def digs(y):
    hn = max((len(str(abs(a.numerator))) for a in y), default=1)
    hd = max((len(str(a.denominator)) for a in y), default=1)
    return hn, hd
hnr, hdr = digs(ye['y_rat']); hns, hds = digs(ye['y_sq'])
cert = {'sector': list(sec), 'rank': len(S), 'ncols': d2['ncols'], 'nrows': d2['nrows'],
        'p1': d2['p'], 'p2': d2.get('p2'), 'dixon_steps': ye['steps'], 'Ls': {k: v for k, v in ye['Ls'].items()},
        'heights': {'rat_num_digits': hnr, 'rat_den_digits': hdr, 'sq_num_digits': hns, 'sq_den_digits': hds},
        'family_counts_used': famcount, 'rows': cert_rows,
        'note': 'sum_j y_j * tuple_j == fold_sector(D) exactly (both components); tuple_j = (1/lam)*fold of src row'}
with open(os.path.join(HERE,f'certificate_{secname}.json'),'w') as f:
    json.dump(cert, f)
log(f'step_04 {secname}: certificate_{secname}.json written rows={len(cert_rows)} bytes={os.path.getsize(os.path.join(HERE,f"certificate_{secname}.json"))} ({time.monotonic()-t0:.0f}s)')
