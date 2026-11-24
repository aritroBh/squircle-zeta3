"""W5G3R: own exact verifier of (I5) and (I4) for the stored certificates, written from FORMAT.md only.
Exact over Z: per target all rationals are brought to one common denominator D (lcm of all denominators),
both sides are multiplied by D and compared as Python integers on every word.
Usage: verify.py [name ...]   (no names: all 74 targets; finished targets in verify.log are skipped)"""
from common import *
from math import lcm
LOG = 'verify.log'
FAMS = ('std', 'lift1', 'lift2old', 'lift4old')

class Rows:
    """row data read as data from the pickles named in FORMAT.md"""
    def __init__(self):
        W5 = WORDS(5); I5 = INDEX(5); W4 = WORDS(4); I4 = INDEX(4)
        self.lit = {}
        for fam in FAMS:
            R = load_data(T + f'/W5C1_certificate_scoping/rels_{fam}.pkl')
            rows = []
            for lab, row in R:
                assert lab[0] == fam
                rows.append(tuple((int(k), int(c)) for k, c in row))
                assert all(0 <= k < len(W5) for k, c in row)
            self.lit[fam] = rows
        V = load_data(T + '/W5C1_certificate_scoping/validated_w5.pkl')
        self.framew5 = [tuple((I5[tuple(w)], int(c)) for w, c in row) for row in V['rows']]   # KeyError = non-admissible
        self.framew5_words = [tuple((tuple(w), int(c)) for w, c in row) for row in V['rows']]
        self.frame234 = load_data(T + '/W5C1_certificate_scoping/rels_frame234.pkl')
        TT = load_data(T + '/A3h_weight4_target/w4_systems_T.pkl')['T0']
        S = load_data(T + '/A3h5_exact_certificate/stage2_rows_S.pkl')['S']
        self.g4S = []
        for i in S:
            d = []
            for w, (q, c) in TT[i][1]:
                assert q == 0
                d.append((I4[tuple(w)], Fraction(c)))
            self.g4S.append(d)
        self.a3h4 = {n: load_data(T + f'/A3h4_frame_ds/validated_w{n}.pickle')['rows'] for n in (2, 3, 4)}

def a3h4_lift(rows, n, rowidx, word):
    """G = A3h4 row (weight n) shuffled with word (weight 4-n), as dict word index -> int (not normalised)"""
    I4 = INDEX(4)
    d = shuffle_vec({tuple(w): int(c) for w, c in rows.a3h4[n][rowidx]}, {tuple(word): 1})
    return {I4[w]: c for w, c in d.items()}

def verify(name, R, corrupt=None, v_override=None, quiet=False):
    """returns dict with I5/I4 mismatch counts. corrupt: ('lit'|'W'|'frame', position) adds 1 to one coefficient."""
    t0 = time.time(); res = {'name': name}
    W5 = WORDS(5); I5 = INDEX(5); W4 = WORDS(4); I4 = INDEX(4); P7 = PERM(5, 7)
    core = read_jsonl_xz(f'{STORE}/{name}.core.jsonl.xz')
    lit = read_jsonl_xz(f'{STORE}/{name}.lit.jsonl.xz')
    t_read = time.time() - t0
    sec = {}
    for o in core:
        assert o['section'] not in sec, ('duplicate section', o['section']); sec[o['section']] = o
    assert set(sec) == {'meta', 'frame', 'W', 't'}, set(sec)
    meta = sec['meta']; assert meta['name'] == name and meta['part'] == 'core'
    u, i = meta['u'], meta['i']; assert name == f'u{u}.f{i}' and u in (3, 5)
    lmeta = lit[0]; assert lmeta['section'] == 'meta' and lmeta['part'] == 'lit' and lmeta['name'] == name
    litrows = {}
    for o in lit[1:]:
        assert o['section'] == 'rows' and o['family'] in FAMS and o['family'] not in litrows
        litrows[o['family']] = o['terms']
    # --- parse every rational; collect denominators
    frame = [(t[:-1], rat(t[-1])) for t in sec['frame']['terms']]
    Wt = [(tuple(w), rat(c)) for w, c in sec['W']['terms']]
    lt = {fam: [(int(k), int(g), rat(c)) for k, g, c in terms] for fam, terms in litrows.items()}
    tt = [(t[:-1], rat(t[-1])) for t in sec['t']['terms']]
    # counts vs meta
    res['n_frame'] = len(frame); res['n_rows'] = {f: len(v) for f, v in lt.items()}; res['W_support'] = len(Wt)
    assert len(frame) == meta['n_frame'] and len(Wt) == meta['W_support']
    assert res['n_rows'] == {f: c for f, c in meta['n_rows'].items() if c}, (res['n_rows'], meta['n_rows'])
    if corrupt:
        kind, pos = corrupt
        if kind == 'lit':
            fam = 'std'; k, g, (p, q) = lt[fam][pos]; lt[fam][pos] = (k, g, (p + q, q))
        elif kind == 'W':
            w, (p, q) = Wt[pos]; Wt[pos] = (w, (p + q, q))
        elif kind == 'frame':
            lab, (p, q) = frame[pos]; frame[pos] = (lab, (p + q, q))
    # ---------------- (I5) exact over Z
    D = 1
    for _, (p, q) in frame: D = lcm(D, q)
    for _, (p, q) in Wt: D = lcm(D, q)
    for terms in lt.values():
        for _, _, (p, q) in terms: D = lcm(D, q)
    res['D5_digits'] = len(str(D))
    acc = [0] * len(W5)
    gs = {}
    for fam, terms in lt.items():
        rows = R.lit[fam]
        seen = set()
        for k, g, (p, q) in terms:
            assert g in (1, 7), g
            assert (k, g) not in seen; seen.add((k, g))
            gs[(fam, g)] = gs.get((fam, g), 0) + 1
            C = p * (D // q)
            if g == 1:
                for wi, c in rows[k]: acc[wi] += C * c
            else:
                for wi, c in rows[k]: acc[P7[wi]] += C * c
    res['g_counts'] = {f'{a}:{b}': c for (a, b), c in sorted(gs.items())}
    t_lit = time.time() - t0
    nf = {'framew5': 0, 'frame234': 0}
    for lab, (p, q) in frame:
        C = p * (D // q)
        if lab[0] == 'framew5':
            _, j = lab
            for wi, c in R.framew5[j]: acc[wi] += C * c
        elif lab[0] == 'frame234':
            _, idx, n, rowidx, word = lab
            flab, frow = R.frame234[idx]
            assert flab == ('frame', n, rowidx, tuple(word)), (flab, lab)
            for wi, c in frow: acc[wi] += C * c
        else: raise ValueError(lab)
        nf[lab[0]] += 1
    res['frame_kinds'] = nf
    Wvec = {}
    for w, (p, q) in Wt:
        assert len(w) == 4 and w in I4 and w not in Wvec, w      # admissible weight-4 word, no duplicates
        Wvec[w] = p * (D // q)
    for w4, c4 in Wvec.items():
        for e, ce in E.items():
            for w, m in shuffle(e, w4).items(): acc[I5[w]] += ce * c4 * m
    # target v = act(framew5 row i, u), times D
    tgt = [0] * len(W5)
    vrow = R.framew5_words[i] if v_override is None else v_override
    for w, c in vrow: tgt[I5[act(w, u)]] += D * c
    mism5 = sum(1 for a, b in zip(acc, tgt) if a != b)
    res['I5_mismatches'] = mism5
    res['I5_target_support'] = sum(1 for b in tgt if b)
    t5 = time.time() - t0
    # ---------------- (I4) exact over Z
    D4 = 1
    for _, (p, q) in Wt: D4 = lcm(D4, q)
    for _, (p, q) in tt: D4 = lcm(D4, q)
    res['D4_digits'] = len(str(D4))
    acc4 = [0] * len(W4); seen = set(); nk = {'A3h5_S': 0, 'A3h4_frame': 0}
    for lab, (p, q) in tt:
        C = p * (D4 // q)
        key = tuple(x if not isinstance(x, list) else tuple(x) for x in lab)
        assert key not in seen; seen.add(key)
        if lab[0] == 'A3h5_S':
            for wi, c in R.g4S[lab[1]]:
                cc = C * c; assert cc.denominator == 1; acc4[wi] += int(cc)
        elif lab[0] == 'A3h4_frame':
            _, n, rowidx, word = lab
            assert n in (2, 3, 4) and len(word) == 4 - n
            for wi, c in a3h4_lift(R, n, rowidx, tuple(word)).items(): acc4[wi] += C * c
        else: raise ValueError(lab)
        nk[lab[0]] += 1
    res['t_kinds'] = nk
    tgt4 = [0] * len(W4)
    for w, (p, q) in Wt: tgt4[I4[w]] += p * (D4 // q)
    res['I4_mismatches'] = sum(1 for a, b in zip(acc4, tgt4) if a != b)
    res['W_nonzero'] = sum(1 for b in tgt4 if b)
    res['time'] = dict(read=round(t_read, 1), lit=round(t_lit - t_read, 1), I5=round(t5, 1), total=round(time.time() - t0, 1))
    return res

if __name__ == '__main__':
    names = sys.argv[1:] or TARGETS
    done = set()
    if os.path.exists(os.path.join(HERE, LOG)):
        for line in open(os.path.join(HERE, LOG)):
            if ' RESULT ' in line: done.add(json.loads(line.split(' RESULT ', 1)[1])['name'])
    todo = [n for n in names if n not in done]
    log(f'start pid {os.getpid()}: {len(todo)} to do, {len(done)} already in log', LOG)
    R = Rows(); MAN = manifest()
    for name in todo:
        # manifest check of the two files read
        ok = all(sha256_file(f'{STORE}/{name}.{p}.jsonl.xz') == MAN[f'{name}.{p}.jsonl.xz'][1] for p in ('core', 'lit'))
        r = verify(name, R); r['sha256_ok'] = ok
        log('RESULT ' + json.dumps(r), LOG)
    log('finished', LOG)
