"""k4h_driver.py -- K4h parallel chamber driver.

Parallelises k4c_driver.chamber_S over the OUTER phi-grid nodes: the outer
node list is built exactly as chamber_S builds it (same splits, same
tanh_sinh_nodes calls, same WMIN filter), cut into contiguous chunks, and
each chunk is accumulated independently; chunk partials are summed in chunk
index order. The per-node inner quadrature (inner_sum) is a VERBATIM copy of
the inner block of k4c_driver.chamber_S, so the only difference vs the
single-process run is summation order across chunk boundaries (roundoff
level).

Multiprocessing layout: init_worker(cfg) runs once per worker (sets dps,
builds model + outer node list); do_chunk(i) processes nodes
[i*chunk, (i+1)*chunk); do_nodes(idxs) processes an arbitrary index list
(used by the pilot). Checkpoint files are JSON with acc as decimal strings
(str(mpf) round-trips exactly -- tested in step_03).
"""
import sys, os, time
K4E = "<project>/research/08_theory/K4e_tetrakis_d3"
if K4E not in sys.path:
    sys.path.insert(0, K4E)
import mpmath as mp
from mpmath import mpf
from k4_engine import FastModel
from k4_quad import tanh_sinh_nodes
from k4c_driver import outer_splits

PI = mp.pi

def build_outer_nodes(normals, HO, TO, WMIN):
    """outer node list, same construction/order as k4c_driver.chamber_S."""
    splits = [mpf(0)] + outer_splits(normals) + [PI/4]
    nodes = []
    for si in range(len(splits)-1):
        phi0, phi1 = splits[si], splits[si+1]
        for (ph, w) in tanh_sinh_nodes(phi0, phi1, HO, TO, WMIN):
            nodes.append((ph, w))
    return nodes

def inner_sum(M, ph, w, K, HI, TI, HI2, WMIN, LMIN):
    """VERBATIM inner block of k4c_driver.chamber_S for one outer node (ph, w).
    Returns (acc[0..K], stats)."""
    m = M.m
    acc = [mpf(0)]*(K+1)
    stats = {'zero_len': 0, 'nbp': 0, 'nsing': 0, 'ninner': 0}
    cph, sph = mp.cos(ph), mp.sin(ph)
    if sph <= 0:
        return acc, stats
    thlo = mp.atan(1/sph)
    if thlo >= PI/2 - LMIN:
        stats['zero_len'] += 1
        return acc, stats
    bps = []
    for j in range(m):
        nx, ny, nz = M.n[j]
        den = nx*cph + ny*sph
        if den == 0: continue
        t = -nz/den
        if t > 0:
            th = mp.atan(t)
            if th > thlo + LMIN and th < PI/2 - LMIN:
                bps.append(th)
    bps.sort()
    merged = []
    for b in bps:
        if merged and b - merged[-1] <= LMIN: continue
        merged.append(b)
    bounds = [thlo] + merged + [PI/2]
    stats['nbp'] += len(merged)
    for i in range(len(bounds)-1):
        lo, hi = bounds[i], bounds[i+1]
        if hi - lo <= LMIN:
            stats['zero_len'] += 1
            continue
        huse = HI2 if (len(merged) > 0) else HI
        if len(merged) > 0: stats['nsing'] += 1
        mid = (lo+hi)/2
        sth_, cth_ = mp.sin(mid), mp.cos(mid)
        ux_, uy_, uz_ = sth_*cph, sth_*sph, cth_
        sig = tuple(1 if M.n[j][0]*ux_+M.n[j][1]*uy_+M.n[j][2]*uz_ > 0 else -1 for j in range(m))
        def g(th):
            sth, cth = mp.sin(th), mp.cos(th)
            ux, uy, uz = sth*cph, sth*sph, cth
            Fs = M.F_series(ux, uy, uz, sig)
            return [Fs[k]*sth for k in range(K+1)]
        nl = tanh_sinh_nodes(lo, hi, huse, TI, WMIN)
        stats['ninner'] += len(nl)
        for (x, w2) in nl:
            gx = g(x)
            for k in range(K+1):
                acc[k] += w*w2*gx[k]
    return acc, stats

# ---------------- worker state (one per process) ----------------
_CFG = None
_NODES = None
_MODEL = None

def init_worker(cfg):
    global _CFG, _NODES, _MODEL
    mp.mp.dps = int(cfg['dps']) + int(cfg.get('guard', 0))
    normals = [[mpf(x) for x in row] for row in cfg['normals']]
    c = dict(cfg)
    for key in ('WMIN', 'LMIN', 'HO', 'HI', 'HI2', 'TO', 'TI'):
        c[key] = mpf(cfg[key])
    _CFG = c
    _MODEL = FastModel(normals, int(cfg['K']))
    _NODES = build_outer_nodes(normals, c['HO'], c['TO'], c['WMIN'])

def _accum(idxs):
    K = int(_CFG['K'])
    acc = [mpf(0)]*(K+1)
    stats = {'zero_len': 0, 'nbp': 0, 'nsing': 0, 'ninner': 0}
    for i in idxs:
        ph, w = _NODES[i]
        a, st = inner_sum(_MODEL, ph, w, K, _CFG['HI'], _CFG['TI'], _CFG['HI2'],
                          _CFG['WMIN'], _CFG['LMIN'])
        for k in range(K+1):
            acc[k] += a[k]
        for key in stats:
            stats[key] += st[key]
    return acc, stats

def do_chunk(idx):
    t0 = time.time()
    CS = int(_CFG['chunk'])
    lo = idx*CS
    hi = min(len(_NODES), lo+CS)
    acc, stats = _accum(range(lo, hi))
    return {'idx': int(idx), 'n': hi-lo,
            'acc': [mp.nstr(v, int(_CFG['dps'])+int(_CFG.get('guard',0))+10) for v in acc],
            'stats': stats, 'el': time.time()-t0}

def do_nodes(idxs):
    """pilot worker: arbitrary index list."""
    t0 = time.time()
    acc, stats = _accum(list(idxs))
    return {'idx': -1, 'n': len(idxs),
            'acc': [mp.nstr(v, int(_CFG['dps'])+int(_CFG.get('guard',0))+10) for v in acc],
            'stats': stats, 'el': time.time()-t0}

# ---------------- combine ----------------
def combine(chunk_results, K):
    """sum chunk partials in chunk index order; returns (acc, stats)."""
    acc = [mpf(0)]*(K+1)
    stats = {'zero_len': 0, 'nbp': 0, 'nsing': 0, 'ninner': 0}
    for r in sorted(chunk_results, key=lambda r: r['idx']):
        for k in range(K+1):
            acc[k] += mpf(r['acc'][k])
        for key in stats:
            stats[key] += r['stats'].get(key, 0)
    return acc, stats

def chamber_S_par(parent_nodes_count_only=False):
    raise NotImplementedError  # driver is used via job scripts (see step_04/jobA_h768)


# ---------------- chunked run with checkpoints (shared by validation + production) ----------------
import glob, json
import multiprocessing as M

def run_job(cfg, ckdir, workers, log=print, tag="job"):
    """Run all chunks of the outer grid with a spawn Pool, checkpointing each
    chunk partial to ckdir/chunk_XXXX.json (atomic via tmp+replace). Resumes:
    chunks already present are skipped. Returns (acc, stats, N, nch)."""
    mp.mp.dps = int(cfg['dps']) + int(cfg.get('guard', 0))
    normals = [[mpf(x) for x in row] for row in cfg['normals']]
    nodes = build_outer_nodes(normals, mpf(cfg['HO']), mpf(cfg['TO']), mpf(cfg['WMIN']))
    N = len(nodes)
    CS = int(cfg['chunk'])
    nch = (N + CS - 1)//CS
    os.makedirs(ckdir, exist_ok=True)
    done = {}
    for f in sorted(glob.glob(os.path.join(ckdir, "chunk_*.json"))):
        try:
            r = json.load(open(f))
            done[int(r['idx'])] = r
        except Exception:
            log("[%s] WARNING: corrupt checkpoint %s ignored (will recompute)" % (tag, f))
    c0 = len(done)
    pending = [i for i in range(nch) if i not in done]
    log("[%s] outer nodes=%d chunks=%d (size=%d) resumedone=%d pending=%d"
        % (tag, N, nch, CS, c0, len(pending)))
    t0 = time.time()
    if pending:
        ctx = M.get_context("spawn")
        with ctx.Pool(int(workers), initializer=init_worker, initargs=(cfg,)) as pool:
            for r in pool.imap_unordered(do_chunk, pending):
                tmp = os.path.join(ckdir, "chunk_%04d.json.tmp" % r['idx'])
                json.dump(r, open(tmp, "w"))
                os.replace(tmp, os.path.join(ckdir, "chunk_%04d.json" % r['idx']))
                done[int(r['idx'])] = r
                el = time.time() - t0
                nthis = len(done) - c0
                rate = nthis/el if el > 0 else 0.0
                eta = (len(pending) - nthis)/rate if rate > 0 else float('inf')
                log("[%s] chunk %d/%d done n=%d el=%.0fs | progress %d/%d elapsed=%.0fs eta=%.0fs"
                    % (tag, r['idx'], nch-1, r['n'], r['el'], len(done), nch, el, eta))
    assert len(done) == nch, "missing chunks: %d/%d" % (len(done), nch)
    K = int(cfg['K'])
    acc, stats = combine(list(done.values()), K)
    return acc, stats, N, nch
