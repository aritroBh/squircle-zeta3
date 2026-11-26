"""r04 (referee, own code): the closed form of M_5(alpha) (term lists of the claim, with real zeta values),
evaluated by an own HPL evaluator (fixed-point integer Taylor sums of H_w(+-c) from the definition), against
M_5 from the definition (r02, tanh-sinh on the u-Taylor coefficient of the integrand) at two working precisions.
usage: python r04_numeric_compare.py DPS1 DPS2 alpha ...   Does not import m5fit/g8_allN/g6_common."""
import sys, json, re
import mpmath as mp
import r02_definition as R

PBITS = 400

def parse(fname):
    out = []
    for nm, v in json.load(open(fname)):
        q = mp.mpf(mp.fraction(*[int(t) for t in (v.split('/') + ['1'])[:2]]))
        if nm == '1': out.append(('const', 0, None, q)); continue
        m = re.match(r'c\^(-?\d+)$', nm)
        if m: out.append(('const', int(m.group(1)), None, q)); continue
        m = re.match(r'c\^(-?\d+) (He?)\[(.*)\]$', nm)
        out.append((m.group(2), int(m.group(1)), tuple(int(t) for t in m.group(3).split(',')), q))
    return out

def hpl_values(words, x):
    """{w: H_w(x)} for all words, |x|<1, fixed point 2^PBITS, depth-first over suffixes."""
    one = 1 << PBITS
    X = int(mp.floor(x*one))
    N = int(PBITS*0.6931471805599453/(-mp.log(abs(x)))) + 200
    allw = set()
    for w in words:
        for i in range(len(w)): allw.add(w[i:])
    children = {}
    for w in allw:
        if len(w) > 1: children.setdefault(w[1:], []).append(w)
    res = {}
    def visit(w, b):
        res[w] = mp.mpf(sum(b))/one
        for ch in children.get(w, []):
            a = ch[0]
            nb = [0]*(N+1)
            if a == 0:
                for m in range(1, N+1): nb[m] = b[m]//m if b[m] >= 0 else -((-b[m])//m)
            else:
                S = 0; aX = a*X
                for m in range(1, N+1):
                    S = (aX*S >> PBITS) + b[m-1]
                    nb[m] = (X*S >> PBITS)//m
            visit(ch, nb)
    for a in (1, -1):
        if (a,) in allw:
            b = [0]*(N+1); p = one
            for m in range(1, N+1):
                p = (p*X) >> PBITS             # x^m
                b[m] = (p if a == 1 else (p if m % 2 else -p))//m   # a^(m+1) x^m / m ... H_1: x^m/m ; H_-1: (-1)^(m+1) x^m/m
            visit((a,), b)
    return res

def eval_terms(terms, c, Hp, Hm):
    tot = mp.mpf(0)
    for kind, e, w, q in terms:
        if kind == 'const': tot += q*c**e
        elif kind == 'He': tot += q*(c**e*Hp[w] + (-c)**e*Hm[w])
        else: tot += q*c**e*Hp[w]
    return tot

if __name__ == '__main__':
    d1, d2 = int(sys.argv[1]), int(sys.argv[2])
    D = "claim_copy/"
    names = ['1', 'z2', 'z2^2', 'z3', 'z4', 'z5', 'z2*z3']
    mp.mp.dps = d2 + 10
    parts = {nm: parse(D + f"hpl_fit_k5_{nm}.json") for nm in names}
    words = {t[2] for nm in names for t in parts[nm] if t[0] != 'const'}
    for a in sys.argv[3:]:
        mp.mp.dps = d2 + 10
        al = mp.mpmathify(eval(a, {'pi': mp.pi}))
        c = mp.cos(al/2)
        Hp = hpl_values(words, c); Hm = hpl_values(words, -c)
        z2, z3, z4, z5 = mp.zeta(2), mp.zeta(3), mp.zeta(4), mp.zeta(5)
        V = {nm: eval_terms(parts[nm], c, Hp, Hm) for nm in names}
        closed = V['1'] + z2*V['z2'] + z3*V['z3'] + z4*V['z4'] + z2**2*V['z2^2'] + z5*V['z5'] + z2*z3*V['z2*z3']
        mp.mp.dps = d1; v1 = R.Mk(al)
        mp.mp.dps = d2; v2 = R.Mk(al)
        mp.mp.dps = d2 + 10
        print(f"alpha={a} c={mp.nstr(c, 12)}: definition(dps {d2}) = {mp.nstr(v2, d2-3)}", flush=True)
        print(f"    |def(dps {d1}) - def(dps {d2})| = {mp.nstr(abs(v1-v2), 3)};  |closed form - def(dps {d2})| = {mp.nstr(abs(closed-v2), 3)};  relative {mp.nstr(abs(closed-v2)/abs(v2), 3)}", flush=True)
