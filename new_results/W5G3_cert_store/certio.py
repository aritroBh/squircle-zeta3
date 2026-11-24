"""W5G3 certificate store writer. One gzip file per target, JSON Lines, one JSON object per line (see FORMAT.md).
t04 creates the file (sections meta, frame, rows (one line per family), W); t05 appends the section t as a further
gzip member (concatenated gzip members read back as one stream)."""
import gzip, lzma, json, os
from fractions import Fraction as F
STORE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'store')
def fstr(x):
    x = F(x); return f'{x.numerator}/{x.denominator}'
def jl(obj): return (json.dumps(obj, separators=(',', ':')) + '\n').encode()
def target_name(i, u): return f'u{u}.f{i}'
def path(i, u, part): return os.path.join(STORE, target_name(i, u) + f'.{part}.jsonl.xz')
def write_lines(p, gen, mode):
    with lzma.open(p, mode, preset=3) as f:      # mode 'ab' appends a further xz stream (read back transparently)
        for b in gen: f.write(b)
def lines_core(i, u, frame, W, extra):
    yield jl(dict(section='meta', name=target_name(i, u), part='core', i=i, u=u, format=2, **extra))
    yield jl(dict(section='frame', terms=frame))
    yield jl(dict(section='W', terms=W))
def lines_lit(i, u, rows, extra):
    yield jl(dict(section='meta', name=target_name(i, u), part='lit', i=i, u=u, format=2, **extra))
    for fam in ('std', 'lift1', 'lift2old', 'lift4old'):
        if fam in rows: yield jl(dict(section='rows', family=fam, terms=rows[fam]))
def lines_w5(i, u, frame, rows, W, extra):
    """frame: list of lists; rows: {family: list of [idx, g, 'p/q']}; W: list of [word, 'p/q']"""
    yield jl(dict(section='meta', name=target_name(i, u), i=i, u=u, format=1, **extra))
    yield jl(dict(section='frame', terms=frame))
    for fam in ('std', 'lift1', 'lift2old', 'lift4old'):
        if fam in rows:
            yield jl(dict(section='rows', family=fam, terms=rows[fam]))
    yield jl(dict(section='W', terms=W))

def verify_pair(i, u):
    """Full integrity check of the two store files of target (i,u): complete xz decode (CRC64 checked by xz), every JSON line parses,
    sections present, term counts equal the counts recorded in meta. Returns (ok, info, frame_terms, W, sha256_core, sha256_lit)."""
    import hashlib
    try:
        res = {}
        for part in ('core', 'lit'):
            p = path(i, u, part)
            if not os.path.exists(p) or os.path.getsize(p) == 0: return False, f'{part} missing/empty', None, None, None, None
            h = hashlib.sha256(); secs = []
            with lzma.open(p, 'rb') as f:
                for line in f:
                    h.update(line); secs.append(json.loads(line))
            res[part] = (secs, h.hexdigest())
        core, lit = res['core'][0], res['lit'][0]
        m = core[0]; ml = lit[0]
        assert m['section'] == 'meta' and m['part'] == 'core' and m['name'] == target_name(i, u), 'core meta'
        assert ml['section'] == 'meta' and ml['part'] == 'lit' and ml['name'] == target_name(i, u), 'lit meta'
        fr = [s for s in core if s['section'] == 'frame']; Wl = [s for s in core if s['section'] == 'W']
        assert len(fr) == 1 and len(Wl) == 1 and len(fr[0]['terms']) == m['n_frame'] and len(Wl[0]['terms']) == m['W_support'], 'core counts'
        rws = {s['family']: len(s['terms']) for s in lit if s['section'] == 'rows'}
        assert rws == ml['n_rows'] == m['n_rows'], 'lit counts'
        W = {tuple(w): F(*map(int, s.split('/'))) for w, s in Wl[0]['terms']}
        return True, 'ok', fr[0]['terms'], W, res['core'][1], res['lit'][1]
    except Exception as e:
        return False, f'{type(e).__name__}: {e}', None, None, None, None
