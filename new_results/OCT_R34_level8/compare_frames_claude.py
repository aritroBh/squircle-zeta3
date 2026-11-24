"""Is the new weight-3 direction from the octagon R3/R4 identities the same as the one from A3h4's frame rows?
Ranks mod two primes of base+Kummer (544) with: frame_w3 rows; DIFF_Lam; both; and + R1..R5."""
import sys, pickle
sys.path.insert(0, '<project>/paper/squircle_note/scripts')
exec(open('OCT_R34_level8/membership_ext_claude.py').read().split("F5 = pickle.load")[0])
from qsqrt2 import QS2
FR = pickle.load(open('A3h4_frame_ds/validated_w3.pickle', 'rb'))['rows']
F5 = pickle.load(open('A3f4_glanois_basis/five_relations.pickle', 'rb'))
for p in (1000000033, 1000000241):
    s = int(sqrt_mod(2, p))
    def mc(a):
        if isinstance(a, QS2): return (modp(a.p, p) + s * modp(a.q, p)) % p
        return modp(a, p)
    def vec(pairs):
        v = [0] * n
        for w, a in pairs: v[idx[tuple(w)]] = (v[idx[tuple(w)]] + mc(a)) % p
        return v
    base = [vec(r.items()) for r in Wrows] + [vec([(w, QS2(x, y)) for w, (x, y) in r.items()]) for r in K]
    frames = [vec(r) for r in FR]
    dl = vec(DIFFS['Lam'].items()); ds = vec(DIFFS['swap'].items())
    Rv = [vec(zip(F5['freewords'], rel)) for rel in F5['relations']]
    rk = lambda rows: flint.nmod_mat(len(rows), n, [x for r in rows for x in r], p).rank()
    print(f"p={p}: base+K {rk(base)} | +frames {rk(base + frames)} | +DIFF_Lam {rk(base + [dl])} | "
          f"+frames+DIFF_Lam {rk(base + frames + [dl])} | +frames+both DIFFs {rk(base + frames + [dl, ds])} | "
          f"+frames+DIFFs+R1..R5 {rk(base + frames + [dl, ds] + Rv)}", flush=True)
