"""W4RX part (a): the exact conversion of W_4 and of every term of R, validated numerically.
Independent code (w4rx_core only).  Reads as DATA: paper/ieee/main.tex is transcribed by hand
below; A3_W4_proof/gpl_expression.tsv and A3c_claude_elementary/Q/vectors.json and
data/original_rows.json are read only to COMPARE."""
import sys, json, re, time, pickle
sys.dont_write_bytecode = True
from fractions import Fraction as F
import mpmath as mp
from w4rx_core import *

T0 = time.time()
def log(s): print(f'[{time.time()-T0:6.1f}s] {s}', flush=True)
ROOT = '<project>/research/08_theory'

# ------------------------------------------------------------------ 1. GPL expression of W_4
# l1 = log((v-b)(v+a)/(4v)) = G1(b;v)+G1(-a;v)-G1(0;v)-lambda   since (1-b)(1+a)=2
# l2 = log((a-v)(v+b)/(4v)) = G1(a;v)+G1(-b;v)-G1(0;v)-lambda   since (a-1)(1+b)=2
# Keys (word, k): G_1(word; v) * lambda^k, G_1(p;v) = log((v-p)/(1-p)).
assert (1 - B) * (1 + A) == K(2) and (A - 1) * (1 + B) == K(2)
L1 = {((B,), 0): K(1), ((-A,), 0): K(1), ((K(0),), 0): K(-1), ((), 1): K(-1)}
L2 = {((A,), 0): K(1), ((-B,), 0): K(1), ((K(0),), 0): K(-1), ((), 1): K(-1)}
def xshuffle(u, v):
    if not u: return [v]
    if not v: return [u]
    return [(u[0],) + w for w in xshuffle(u[1:], v)] + [(v[0],) + w for w in xshuffle(u, v[1:])]
def pmul(*fs):
    out = {((), 0): K(1)}
    for f in fs:
        new = {}
        for (u, k), c in out.items():
            for (v, j), d in f.items():
                for w in xshuffle(u, v):
                    new[(w, k + j)] = new.get((w, k + j), K(0)) + c * d
        out = {key: c for key, c in new.items() if c != K(0)}
    return out
# W4 = int_b^a [ l1^3 ((8-7s)/(48 v) - 1/(3(v+a))) + (s/16) l1^2 l2 / v ] dv,  s = sqrt2
GPL = {}
for z, fs, c in [(K(0), (L1, L1, L1), K(F(8, 48), F(-7, 48))), (-A, (L1, L1, L1), K(F(-1, 3))),
                 (K(0), (L1, L1, L2), K(0, F(1, 16)))]:
    for (w, k), d in pmul(*fs).items():
        GPL[((z,) + w, k)] = GPL.get(((z,) + w, k), K(0)) + c * d
GPL = {key: c for key, c in GPL.items() if c != K(0)}
log(f'own GPL expression: {len(GPL)} terms  [W4 = sum c lambda^k (G_1(w;a)-G_1(w;b))]')

# compare with the notes' exported expression (data)
def parse_pt(s):
    s = s.strip()
    if s == '0': return K(0)
    m = re.fullmatch(r'\((-?\d+)([+-]\d+)\*sqrt2\)', s); assert m, s
    return K(int(m.group(1)), int(m.group(2)))
tsv = {}
for line in open(ROOT + '/A3_W4_proof/gpl_expression.tsv'):
    if line.startswith('#') or line.startswith('coefficient'): continue
    p, q, k, word = line.rstrip('\n').split('\t')
    pts = tuple(parse_pt(x) for x in re.findall(r'\(-?\d+[+-]\d+\*sqrt2\)|\b0\b', word))
    tsv[(pts, int(k))] = K(F(p), F(q))
same = set(tsv) == set(GPL) and all(tsv[key] == GPL[key] for key in GPL)
log(f'GPL expression equals A3_W4_proof/gpl_expression.tsv exactly: {same} ({len(tsv)} rows)')
assert same

# ------------------------------------------------------------------ 2. exact conversion
LOG2K = [{(): F(1)}]
for _ in range(3): LOG2K.append(vmul(LOG2K[-1], LOG2))
def to_K(vec): return {w: K.c(c) for w, c in vec.items()}
W = {}
conv_cache = {}
for (w, k), c in GPL.items():
    forms = [(p, None) for p in w]                     # dlog(x-p)
    term = vadd((convert_real_segment(forms, K(1), A), 1), (convert_real_segment(forms, K(1), B), -1))
    term = vmul(term, LOG2K[k])
    W = vadd((W, 1), (to_K(term), c))
log(f'W vector: support {len(W)}')

# R, transcribed by hand from paper/ieee/main.tex eq:R (T3)
s = K(0, 1); h = F(1, 2)
LI = {}
ATOMS = {}
def li(n, z):
    key = ('Li', n, repr(z))
    if key not in ATOMS:
        frames = li_frames(z)
        assert frames, ('no frame', n, z)
        vecs = [li_vector(n, z, fr) for fr in frames]
        agree = all(v == vecs[0] for v in vecs)
        # Choice of frame: the first admissible (p0,p1,pinf,e) in lexicographic order over the
        # point list (inf,0,1,-1,a,b,-a,-b) (the order used by the notes; any admissible frame
        # gives the same number, see the numerical check of all frames below).
        ORDER = [None, K(0), K(1), K(-1), A, B, -A, -B]
        idx = lambda p: [repr(q) for q in ORDER].index(repr(p))
        chosen = min(range(len(frames)), key=lambda i: tuple(idx(p) for p in frames[i]))
        ATOMS[key] = vecs[chosen]
        LI[key] = dict(z=z, n=n, frames=frames, chosen=frames[chosen], all_frames_agree=agree,
                       nvec=len(set(repr(sorted(v.items())) for v in vecs)), vecs=vecs)
        log(f'Li_{n}({z}): {len(frames)} admissible frames, {LI[key]["nvec"]} distinct vectors; chosen frame (p0,p1,pinf,e) = {frames[chosen]}')
    return ATOMS[key]
Lv, lam, pi2, z3 = LOGA(), LOG2, PI2, ZETA3
b, onebysq2 = B, K(0, h)
R_TERMS = [   # (coefficient, list of atom vectors, name)
    (K(1, -573), [li(4, b)], 'Li4(b)'), (K(-1, -573), [li(4, -b)], 'Li4(-b)'),
    (K(0, 56), [li(4, -b*b)], 'Li4(-b^2)'), (K(0, 17), [li(4, onebysq2)], 'Li4(1/sqrt2)'),
    (K(F(-37, 4), F(-17, 16)), [li(4, K(h))], 'Li4(1/2)'),
    (K(10, 86), [li(4, 1 - onebysq2)], 'Li4(1-1/sqrt2)'), (K(-10, 86), [li(4, 2 - s)], 'Li4(2-sqrt2)'),
    (K(40, 7), [Lv, li(3, b)], 'L Li3(b)'), (K(-40, 7), [Lv, li(3, -b)], 'L Li3(-b)'),
    (K(-72), [Lv, li(3, onebysq2)], 'L Li3(1/sqrt2)'), (K(F(63, 8), F(-3, 16)), [Lv, z3], 'L zeta3'),
    (K(6, 13), [lam, li(3, b)], 'lam Li3(b)'), (K(6, -13), [lam, li(3, -b)], 'lam Li3(-b)'),
    (K(0, F(-47, 2)), [lam, li(3, onebysq2)], 'lam Li3(1/sqrt2)'),
    (K(F(-3, 8), F(329, 128)), [lam, z3], 'lam zeta3'),
    (K(30), [li(2, b), li(2, b)], 'Li2(b)^2'),
    (K(30, F(7, 2)), [li(2, b), Lv, Lv], 'Li2(b) L^2'), (K(-24), [li(2, b), lam, Lv], 'Li2(b) lam L'),
    (K(0, F(-19, 8)), [li(2, b), lam, lam], 'Li2(b) lam^2'),
    (K(F(-19, 4), F(1, 24)), [li(2, b), pi2], 'Li2(b) pi^2'),
    (K(F(103, 384), F(-4963, 5760)), [pi2, pi2], 'pi^4'),
    (K(F(67, 24), F(185, 96)), [pi2, Lv, Lv], 'pi^2 L^2'),
    (K(F(-35, 12), F(41, 12)), [pi2, lam, Lv], 'pi^2 lam L'),
    (K(F(5, 12), F(-1105, 384)), [pi2, lam, lam], 'pi^2 lam^2'),
    (K(F(13, 24), F(-25, 24)), [Lv] * 4, 'L^4'), (K(F(-73, 6), F(-29, 3)), [Lv] * 3 + [lam], 'L^3 lam'),
    (K(F(63, 8), F(67, 16)), [Lv, Lv, lam, lam], 'L^2 lam^2'),
    (K(F(41, 24), F(17, 6)), [Lv] + [lam] * 3, 'L lam^3'),
    (K(F(-61, 192), F(137, 192)), [lam] * 4, 'lam^4')]
R = {}
TERMVEC = []
for c, atoms, name in R_TERMS:
    t = {(): F(1)}
    for a in atoms: t = vmul(t, a)
    TERMVEC.append(t)
    R = vadd((R, 1), (to_K(t), c))
D = vadd((W, 1), (R, -1))
def split(vec):
    return [{w: c.p for w, c in vec.items() if c.p}, {w: c.q for w, c in vec.items() if c.q}]
Ws, Rs, Ds = split(W), split(R), split(D)
log(f'supports (rational, sqrt2):  W {tuple(map(len, Ws))}  R {tuple(map(len, Rs))}  D {tuple(map(len, Ds))}')
assert all(len(w) == 4 and admissible(w) for vec in Ds for w in vec)

# ------------------------------------------------------------------ 3. compare with the notes' vectors
notes = json.load(open(ROOT + '/A3c_claude_elementary/Q/vectors.json'))
def notes_vec(lst): return {tuple(w): F(c) for w, c in lst}
cmp = {}
for name, mine in [('W', Ws), ('R', Rs), ('D', Ds)]:
    for j in range(2):
        theirs = notes_vec(notes[name][j]); mineH = to_H(mine[j])
        cmp[f'{name}[{j}]'] = (mineH == theirs)
log(f'exact equality with A3c vectors.json (H notation): {cmp}')
orig = json.load(open('data/original_rows.json'))
cmpD = [to_H(Ds[j]) == notes_vec(orig['D'][j]) for j in range(2)]
log(f'exact equality of D with the verifier target rebuilt by the original code: {cmpD}')

# ------------------------------------------------------------------ 4. numerics
vals = load_values('data/G_values_N900.json', 150)
vals2 = load_values('data/G_values_N1100.json', 150)
dmax = max(abs(vals[w] - vals2[w]) for w in vals)
log(f'G values: N=900 vs N=1100 max difference {mp.nstr(dmax, 5)} over {len(vals)} words')
mp.mp.dps = 150
a_, b_, s_ = 1 + mp.sqrt(2), mp.sqrt(2) - 1, mp.sqrt(2)
Lnum, lamnum = mp.log(a_), mp.log(2)
def l1(v): return mp.log((v - b_) * (v + a_) / (4 * v))
def l2(v): return mp.log((a_ - v) * (v + b_) / (4 * v))
def integrand(v):
    return l1(v)**3 * ((8 - 7*s_) / (48 * v) - 1 / (3 * (v + a_))) + s_ / 16 * l1(v)**2 * l2(v) / v
W4q = mp.quad(integrand, [b_, 1, a_], maxdegree=12)
W4q2 = mp.quad(integrand, [b_, mp.mpf('0.7'), 1, mp.mpf('1.6'), a_], maxdegree=12)
log(f'W4 by quadrature = {mp.nstr(W4q, 40)} (two subdivisions differ by {mp.nstr(abs(W4q-W4q2), 3)})')
direct = {'Li': lambda n, z: mp.polylog(n, z.mpf())}
results = {}
W4conv = num(W, vals)
results['W4: converted vector vs quadrature'] = abs(W4conv - W4q)
log(f'W4 converted vector value = {mp.nstr(W4conv.real, 40)}  |diff vs quad| = {mp.nstr(abs(W4conv - W4q), 5)}')
atom_err = {}
for key, info in LI.items():
    v = num(ATOMS[key], vals); d = mp.polylog(info['n'], info['z'].mpf())
    atom_err[f"Li_{info['n']}({info['z']})"] = abs(v - d)
atom_err['log a'] = abs(num(Lv, vals) - Lnum); atom_err['log 2'] = abs(num(lam, vals) - lamnum)
atom_err['pi^2'] = abs(num(pi2, vals) - mp.pi**2); atom_err['zeta(3)'] = abs(num(z3, vals) - mp.zeta(3))
for key, info in LI.items():
    d = mp.polylog(info['n'], info['z'].mpf())
    e = max(abs(num(v, vals) - d) for v in info['vecs'])
    atom_err[f"Li_{info['n']}({info['z']}) all {len(info['vecs'])} frames"] = e
for k, e in atom_err.items(): log(f'  atom {k:22s} |converted - direct| = {mp.nstr(e, 5)}')
def direct_atom(name):
    return None
num_atoms = {'L': Lnum, 'lam': lamnum, 'pi2': mp.pi**2, 'z3': mp.zeta(3)}
Rdirect = 0; term_err = []
for (c, atoms, name), tv in zip(R_TERMS, TERMVEC):
    # direct value of the product of atoms
    dv = mp.mpf(1)
    for a in atoms:
        for key, info in LI.items():
            if ATOMS[key] is a: dv *= mp.polylog(info['n'], info['z'].mpf()); break
        else:
            dv *= {id(Lv): Lnum, id(lam): lamnum, id(pi2): mp.pi**2, id(z3): mp.zeta(3)}[id(a)]
    Rdirect += c.mpf() * dv
    term_err.append((name, abs(num(tv, vals) - dv)))
for name, e in term_err: log(f'  term {name:20s} |converted - direct| = {mp.nstr(e, 5)}')
results['max term error'] = max(e for n_, e in term_err)
results['R: converted vs direct'] = abs(num(R, vals) - Rdirect)
results['W4(quad) - R(direct)'] = abs(W4q - Rdirect)
results['D: value of exact D vector'] = abs(num(D, vals))
for k, v in results.items(): log(f'{k}: {mp.nstr(v, 5)}')
pickle.dump({'W': W, 'R': R, 'D': D, 'Ds': Ds, 'Ws': Ws, 'Rs': Rs, 'LI': {k: {**{x: y for x, y in v.items() if x != 'vecs'}, 'chosen': tuple(map(repr, v['chosen'])), 'z': repr(v['z']), 'frames': [tuple(map(repr, f)) for f in v['frames']]} for k, v in LI.items()},
             'W4': mp.nstr(W4q, 130)}, open('data/check_a.pickle', 'wb'))
log('saved data/check_a.pickle')
