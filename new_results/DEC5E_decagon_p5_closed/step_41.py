src = '''"""step_40_fits.py - LLL fits of Wt3, Wt4, Wt5 at two precisions over Q(sqrt5).

Basis (task prescription; golden identities already applied):
 ring_k: monomials of degree k in {l2, lf (wt1), z2 (wt2), z3 (wt3)} + z5 at k=5
 extras: L3m = Li3(-1/phi) (wt3); Li4m = Li4(-1/phi), Li4q = Li4(1/phi^2) (wt4);
         Li5m = Li5(-1/phi), Li5q = Li5(1/phi^2) (wt5).
 Wt3 basis: ring3 + L3m                          (8 atoms)
 Wt4 basis: ring4 + Li4m, Li4q, L3m*l2, L3m*lf   (15 atoms)
 Wt5 basis: ring5 + Li5m, Li5q, (Li4m,Li4q)*log, L3m*log^2, z2*L3m  (27 atoms)
Usage: python step_40_fits.py LEVEL DPS1 DPS2
"""
import sys, os, pickle, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from q5_lll import fit, report

def atoms(D):
    mp.mp.dps = D + 30
    phi = (1+mp.sqrt(5))/2
    l2 = mp.log(2); lf = mp.log(phi); z2 = mp.zeta(2); z3 = mp.zeta(3); z5 = mp.zeta(5)
    L3m = mp.polylog(3, -1/phi); L4m = mp.polylog(4, -1/phi); L4q = mp.polylog(4, 1/phi**2)
    L5m = mp.polylog(5, -1/phi); L5q = mp.polylog(5, 1/phi**2)
    A = {'l2': l2, 'lf': lf, 'z2': z2, 'z3': z3, 'z5': z5,
         'L3m': L3m, 'L4m': L4m, 'L4q': L4q, 'L5m': L5m, 'L5q': L5q}
    out = {}
    def add(name, val): out[name] = val
    # ring3
    for a in range(4):
        b = 3-a
        add(f'l2^{a}lf^{b}' if a and b else (f'l2^{a}' if a else f'lf^{b}'), l2**a * lf**b)
    add('z2l2', z2*l2); add('z2lf', z2*lf); add('z3', z3)
    # ring4
    for a in range(5):
        b = 4-a
        add(f'l2^{a}lf^{b}' if a and b else (f'l2^{a}' if a else f'lf^{b}'), l2**a * lf**b)
    for a in range(3):
        b = 2-a
        add(f'z2l2^{a}lf^{b}' if a and b else (f'z2l2^{a}' if a else f'z2lf^{b}'), z2*l2**a*lf**b)
    add('z3l2', z3*l2); add('z3lf', z3*lf); add('z2^2', z2**2)
    # ring5
    for a in range(6):
        b = 5-a
        add(f'l2^{a}lf^{b}' if a and b else (f'l2^{a}' if a else f'lf^{b}'), l2**a * lf**b)
    for a in range(4):
        b = 3-a
        add(f'z2l2^{a}lf^{b}' if a and b else (f'z2l2^{a}' if a else f'z2lf^{b}'), z2*l2**a*lf**b)
    for a in range(3):
        b = 2-a
        add(f'z3l2^{a}lf^{b}' if a and b else (f'z3l2^{a}' if a else f'z3lf^{b}'), z3*l2**a*lf**b)
    add('z2^2l2', z2**2*l2); add('z2^2lf', z2**2*lf); add('z5', z5); add('z2z3', z2*z3)
    # extras
    add('L3m', L3m)
    add('L4m', L4m); add('L4q', L4q); add('L3ml2', L3m*l2); add('L3mlf', L3m*lf)
    add('L5m', L5m); add('L5q', L5q)
    add('L4ml2', L4m*l2); add('L4mlf', L4m*lf); add('L4ql2', L4q*l2); add('L4qlf', L4q*lf)
    add('L3ml2^2', L3m*l2**2); add('L3ml2lf', L3m*l2*lf); add('L3mlf^2', L3m*lf**2)
    add('z2L3m', z2*L3m)
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
    W5 = pickle.load(open('wparts_dec_d520.pkl','rb'))
    W8 = pickle.load(open('wparts_dec_d800.pkl','rb'))
    results = []
    for D, W in [(D1, W5), (D2, W8)]:
        mp.mp.dps = D + 30
        T = mp.mpf(W[f'Wt{level}'])
        A = atoms(D + 30)
        at = [(nm, A[nm]) for nm in basis]
        mp.mp.dps = D + 30
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
        pickle.dump({'level': level, 'rel': r0[2], 'height': r0[0]},
                    open(f'fit_wt{level}.pkl','wb'))

if __name__ == "__main__":
    main()
'''
open('step_40_fits.py','w').write(src)
import subprocess, time
PY = '<project>/.venv/bin/python'
t0=time.time()
r = subprocess.run(['nice', PY, 'step_40_fits.py', '3', '450', '700'], capture_output=True, text=True, cwd='.')
print(r.stdout[-4000:])
if r.stderr: print("STDERR:", r.stderr[-500:])
print("wall %.0fs" % (time.time()-t0))
