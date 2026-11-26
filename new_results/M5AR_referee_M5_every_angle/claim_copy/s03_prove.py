"""s03: PROOF (exact, all n) of the HPL forms of every zeta part of M_5, by identity of canonical Z-sum forms,
exactly as g8_allN.py does for k = 2,3,4 (functions mu_side, hpl_side, low_order imported unchanged from the copy
of g8_allN.py).  For each part z:
   even identity: sum_terms [c^(2n)] term  -  [z]P_{5,n}  == 0  as Q(n)-combination of Z-words Z_w(2n)  (all n >= 1)
   odd identity : sum_terms [c^(2n+1)] term == 0                                                          (all n >= 1)
   low orders   : [c^N] of the term list == 0 for N = -6..1, exact rationals.
Plus: the zeta ladder rungs at k = 5 as identities of Z-forms of the mu sides (independent of any fit), and the
comparison of the independently fitted term lists with the term lists predicted by the ladder."""
import sys, json, re, time
sys.setrecursionlimit(20000)
from fractions import Fraction as Fr
import g8_allN as G
from g8_zsum import K, nn, add_into, lin

def parse(fname_or_list):
    raw = json.load(open(fname_or_list)) if isinstance(fname_or_list, str) else fname_or_list
    out = []
    for nm, v in raw:
        if nm == '1': out.append(('const', (), 0, Fr(v))); continue
        if nm == 'c': out.append(('const', (), 1, Fr(v))); continue
        m = re.match(r'c\^(-?\d+)$', nm)
        if m: out.append(('const', (), int(m.group(1)), Fr(v))); continue
        m = re.match(r'c\^(-?\d+) (He?)\[(.*)\]', nm)
        e = int(m.group(1)); w = tuple(int(t) for t in m.group(3).split(','))
        out.append(('even' if m.group(2) == 'He' else 'plain', w, e, Fr(v)))
    return out

def normalise(terms):
    """canonical dict {('const',e) or (e,w): q} of PLAIN terms q c^e H_w(c) (an He pair is split into its two plain terms)."""
    d = {}
    def put(k_, q): d[k_] = d.get(k_, 0) + q
    for kind, w, e, q in terms:
        if kind == 'const': put(('const', e), q)
        elif kind == 'plain': put((e, w), q)
        else:
            s = 1 if (e + G.nz(w)) % 2 == 0 else -1
            put((e, w), q); put((e, G.neg(w)), s*q)
    return {k: v for k, v in d.items() if v != 0}

def check_part(name, terms, muE):
    t0 = time.time()
    ev_even = G.hpl_side(terms, 0)
    diff = lin((1, ev_even), (-1, muE))
    odd = G.hpl_side(terms, 1)
    lo = G.low_order(terms, lambda N: Fr(0), lo=-6, hi=1)
    from collections import Counter
    cw = Counter(len(w) for kind, w, e, q in terms if kind != 'const')
    es = sorted({e for kind, w, e, q in terms if kind != 'const'})
    nconst = sum(1 for t in terms if t[0] == 'const')
    dens = sorted({q.denominator for *_, q in terms})
    print(f"part {name}: {len(terms)} terms ({nconst} constant); terms by weight {dict(sorted(cw.items()))}; e in {es}; lcm-relevant denominators of q: max {max(dens)}")
    print(f"   HPL side {len(ev_even)} words, mu side {len(muE)} words")
    print(f"   even identity (all n>=1): {len(diff) == 0}")
    print(f"   odd identity  (all n>=1): {len(odd) == 0}")
    print(f"   low orders N=-6..1 exact: {len(lo) == 0}   [{time.time()-t0:.0f}s]", flush=True)
    if diff: print('   residual:', list(diff.items())[:3])
    if lo: print('   low mismatch:', lo)
    return not diff and not odd and not lo

if __name__ == '__main__':
    t0 = time.time()
    P5 = G.mu_side(5); P4 = G.mu_side(4); P3 = G.mu_side(3); P2 = G.mu_side(2)
    print(f"mu sides k=2..5 built in {time.time()-t0:.0f}s; k=5 monomials: " +
          ", ".join(f"{G.zmon_name(zm)}({len(E)} words)" for zm, E in sorted(P5.items())), flush=True)
    Z = lambda **kw: tuple(kw.get(f"z{w}", 0) for w in range(2, 6))
    # ---- content theorem: the monomials that occur
    allowed = {Z(), Z(z2=1), Z(z3=1), Z(z4=1), Z(z2=2), Z(z5=1), Z(z2=1, z3=1)}
    print("content: monomials of P_{5,n} =", sorted(G.zmon_name(z) for z in P5), "; all of weight <= 5 and all 7 present:", set(P5) == allowed)
    # ---- ladder rungs at k=5 on the mu side (exact, all n)
    def eq(A, B): return len(lin((1, A), (-1, B))) == 0
    z2_4 = P4[(1, 0, 0)]; z22_4 = P4[(2, 0, 0)]; z2_3 = P3[(1, 0)]; z2_2 = P2[(1,)]
    r1 = eq(P5[Z(z3=1)], lin((-1, z2_4)))
    r2 = eq(P5[Z(z4=1)], lin((1, z2_3), (1, P5[Z(z2=2)])))
    r3 = eq(P5[Z(z5=1)], lin((-1, z2_2), (-2, z22_4)))
    r4 = eq(P5[Z(z2=1, z3=1)], lin((-2, z22_4)))
    print("ladder on mu sides (Q(n)-Z-forms, all n):")
    print("   [z3]P5 == -[z2]P4 :", r1)
    print("   [z4]P5 == [z2]P3 + [z2^2]P5 :", r2)
    print("   [z5]P5 == -[z2]P2 - 2[z2^2]P4 :", r3)
    print("   [z2z3]P5 == -2[z2^2]P4 :", r4)
    import sympy as sp
    for nm, zm in (('z5', Z(z5=1)), ('z2*z3', Z(z2=1, z3=1)), ('z2^2', Z(z2=2)), ('z4', Z(z4=1))):
        print(f"   [{nm}]P5 =", {w: sp.factor(c.as_expr()) for w, c in P5[zm].items()})
    # ---- the all-n identities, each part separately
    files = {'1': Z(), 'z2': Z(z2=1), 'z2^2': Z(z2=2), 'z3': Z(z3=1), 'z4': Z(z4=1), 'z5': Z(z5=1), 'z2*z3': Z(z2=1, z3=1)}
    ok = {}
    for name in ['z2^2', 'z5', 'z2*z3', 'z4', 'z3', 'z2', '1']:
        ok[name] = check_part(name, parse(f"hpl_fit_k5_{name}.json"), P5[files[name]])
    print("ALL-n IDENTITIES:", ok)
    # ---- function level: independent fits vs ladder predictions (term lists, canonical form)
    T = {nm: normalise(parse(f"hpl_fit_k5_{nm}.json")) for nm in files}
    Psi2 = normalise(parse('g8_terms_k4_partz2.json')); Psi = normalise(parse('g8_terms_k3_partz2.json'))
    # Phi = 1/2 - (1-c^2)/(2c) artanh c ; artanh c = He[1]/2 at odd e  ->  1/2 - c^-1 He[1]/4 + c He[1]/4
    Phi = normalise([('const', (), 0, Fr(1, 2)), ('even', (1,), -1, Fr(-1, 4)), ('even', (1,), 1, Fr(1, 4))])
    # F2 = (-c/16 - 1/(8c) + 3/(16c^3)) artanh c + 1/16 - 3/(16 c^2)
    F2 = normalise([('even', (1,), 1, Fr(-1, 32)), ('even', (1,), -1, Fr(-1, 16)), ('even', (1,), -3, Fr(3, 32)), ('const', (), 0, Fr(1, 16)), ('const', (), -2, Fr(-3, 16))])
    def comb(*pairs):
        d = {}
        for a, D in pairs:
            for k_, v in D.items(): d[k_] = d.get(k_, 0) + a*v
        return {k_: v for k_, v in d.items() if v != 0}
    print("ladder predictions vs independent fits (canonical term lists):")
    print("   [z3]M5 == -Psi_2 :", T['z3'] == comb((-1, Psi2)))
    print("   [z4]M5 == Psi + [z2^2]M5 :", T['z4'] == comb((1, Psi), (1, T['z2^2'])))
    print("   [z5]M5 == -Phi - 2 F_2 :", T['z5'] == comb((-1, Phi), (-2, F2)))
    print("   [z2z3]M5 == -2 F_2 :", T['z2*z3'] == comb((-2, F2)))
    print("THEOREM M5 (all parts, all n):", all(ok.values()) and r1 and r2 and r3 and r4)
