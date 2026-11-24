"""step_40_fits.py - LLL fits of Wt3, Wt4, Wt5 at two precisions over Q(sqrt5)."""
import sys, os, pickle, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from q5_lll import fit, report

def mon(base2, a, basef, b):
    if a == 0 and b == 0: return ''
    parts = []
    if a: parts.append(base2 if a == 1 else f'{base2}^{a}')
    if b: parts.append(basef if b == 1 else f'{basef}^{b}')
    return ''.join(parts)

def atoms(D):
    mp.mp.dps = D + 30
    phi = (1+mp.sqrt(5))/2
    l2 = mp.log(2); lf = mp.log(phi); z2 = mp.zeta(2); z3 = mp.zeta(3); z5 = mp.zeta(5)
    L3m = mp.polylog(3, -1/phi); L4m = mp.polylog(4, -1/phi); L4q = mp.polylog(4, 1/phi**2)
    L5m = mp.polylog(5, -1/phi); L5q = mp.polylog(5, 1/phi**2)
    out = {}
    def logs(k):  # pure log monomials of degree k
        for a in range(k+1):
            out[mon('l2', a, 'lf', k-a)] = l2**a * lf**(k-a)
    for k in (3, 4, 5): logs(k)
    for k in (2, 3):
        for a in range(k+1):
            out['z2' + mon('l2', a, 'lf', k-a)] = z2*l2**a*lf**(k-a)
    for k in (1, 2):
        for a in range(k+1):
            out['z3' + mon('l2', a, 'lf', k-a)] = z3*l2**a*lf**(k-a)
    out['z2l2'] = z2*l2; out['z2lf'] = z2*lf; out['z3'] = z3; out['z5'] = z5
    out['z2^2'] = z2**2; out['z2^2l2'] = z2**2*l2; out['z2^2lf'] = z2**2*lf
    out['z2z3'] = z2*z3
    out['L3m'] = L3m
    out['L4m'] = L4m; out['L4q'] = L4q; out['L3ml2'] = L3m*l2; out['L3mlf'] = L3m*lf
    out['L5m'] = L5m; out['L5q'] = L5q
    out['L4ml2'] = L4m*l2; out['L4mlf'] = L4m*lf; out['L4ql2'] = L4q*l2; out['L4qlf'] = L4q*lf
    out['L3ml2^2'] = L3m*l2**2; out['L3ml2lf'] = L3m*l2*lf; out['L3mlf^2'] = L3m*lf**2
    out['z2L3m'] = z2*L3m
    return out

RING3 = ['l2^3','l2^2lf','l2lf^2','lf^3','z2l2','z2lf','z3']
B3 = RING3 + ['L3m']
RING4 = ['l2^4','l2^3lf','l2^2lf^2','l2lf^3','lf^4','z2l2^2','z2l2lf','z2lf^2','z3l2','z3lf','z2^2']
B4 = RING4 + ['L4m','L4q','L3ml2','L3mlf']
RING5 = ['l2^5','l2^4lf','l2^3lf^2','l2^2lf^3','l2lf^4','lf^5',
         'z2l2^3','z2l2^2lf','z2l2lf^2','z2lf^3',
         'z3l2^2','z3l2lf','z3lf^2','z2^2l2','z2^2lf','z5','z2z3']
B5 = RING5 + ['L5m','L5q','L4ml2','L4mlf','L4ql2','L4qlf','L3ml2^2','L3ml2lf','L3mlf^2','z2L3m']

def main():
    level = int(sys.argv[1]); D1, D2 = int(sys.argv[2]), int(sys.argv[3])
    basis = {3: B3, 4: B4, 5: B5}[level]
    W5 = pickle.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'wparts_dec_d520.pkl'),'rb'))
    W8 = pickle.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'wparts_dec_d800.pkl'),'rb'))
    results = []
    for D, W in [(D1, W5), (D2, W8)]:
        mp.mp.dps = D + 40
        T = mp.mpf(W[f'Wt{level}'])
        A = atoms(D + 10)
        at = [(nm, A[nm]) for nm in basis]
        cands = fit(T, at, cexp=int(D*0.92), nbest=4)
        print(f"--- Wt{level} fit at {D} digits ({len(basis)} atoms) ---")
        report(f"Wt{level}@{D}", cands)
        results.append(cands)
    r0, r1 = results[0][0], results[1][0]
    same = r0[2] == r1[2]
    print(f">>> Wt{level}: identical at {D1} and {D2}: {same}")
    if same:
        print("    residuals: %s / %s ; height %s" % (mp.nstr(mp.mpf(r0[1]),3), mp.nstr(mp.mpf(r1[1]),3), r0[0]))
        rel = {k: (str(v[0]), str(v[1])) for k, v in sorted(r0[2].items())}
        print(json.dumps(rel, indent=1))
        here = os.path.dirname(os.path.abspath(__file__))
        pickle.dump({'level': level, 'rel': r0[2], 'height': r0[0]},
                    open(os.path.join(here, f'fit_wt{level}.pkl'),'wb'))

if __name__ == "__main__":
    main()
