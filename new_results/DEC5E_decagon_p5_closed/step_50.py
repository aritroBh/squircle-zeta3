src = '''"""step_48_rfit.py - fit the symbol-zero remainder R = Wt5 - exact log part.

R basis (21 invisible atoms, all symbol-zero):
 ring-invisible 11: z2l2^3, z2l2^2lf, z2l2lf^2, z2lf^3, z3l2^2, z3l2lf, z3lf^2, z2^2l2, z2^2lf, z5, z2z3
 z2L3m (1)
 K-combos (9, PROVED symbol-zero): K5m, K5q, K4m*l2, K4m*lf, K4q*l2, K4q*lf, K3m*l2^2, K3m*l2lf, K3m*lf^2
Usage: python step_48_rfit.py D1 D2
"""
import sys, os, pickle, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from q5_lll import fit, report

LOGPART = {   # exact coefficients from DEC5D's exact symbol (step_47)
 'l2^5':      ("10", "0"),
 'l2^4lf':    ("125", "-25"),
 'l2^3lf^2':  ("2025/2", "-825/2"),
 'l2^2lf^3':  ("52075/8", "-23525/8"),
 'l2lf^4':    ("449025/16", "-403175/32"),
 'lf^5':      ("15112175/256", "-6759335/256"),
}

def setup(D):
    mp.mp.dps = D + 40
    phi = (1+mp.sqrt(5))/2
    l2 = mp.log(2); lf = mp.log(phi); z2 = mp.zeta(2); z3 = mp.zeta(3); z5 = mp.zeta(5)
    L3m = mp.polylog(3, -1/phi); L4m = mp.polylog(4, -1/phi); L4q = mp.polylog(4, 1/phi**2)
    L5m = mp.polylog(5, -1/phi); L5q = mp.polylog(5, 1/phi**2)
    K5m = L5m + lf**5/120; K5q = L5q - mp.mpf(2)/15*lf**5
    K4m = L4m - lf**4/24;  K4q = L4q + lf**4/3
    K3m = L3m + lf**3/6
    A = {'l2': l2, 'lf': lf, 'z2': z2, 'z3': z3, 'z5': z5}
    atoms = {}
    for a in range(4):
        b = 3-a
        nm = 'z2' + (('l2^%d'%a if a>1 else 'l2'*a) + ('lf^%d'%b if b>1 else 'lf'*b))
        atoms[nm] = z2*l2**a*lf**b
    for a in range(3):
        b = 2-a
        nm = 'z3' + (('l2^%d'%a if a>1 else 'l2'*a) + ('lf^%d'%b if b>1 else 'lf'*b))
        atoms[nm] = z3*l2**a*lf**b
    atoms['z2^2l2'] = z2**2*l2; atoms['z2^2lf'] = z2**2*lf
    atoms['z5'] = z5; atoms['z2z3'] = z2*z3
    atoms['z2L3m'] = z2*L3m
    atoms['K5m'] = K5m; atoms['K5q'] = K5q
    atoms['K4ml2'] = K4m*l2; atoms['K4mlf'] = K4m*lf
    atoms['K4ql2'] = K4q*l2; atoms['K4qlf'] = K4q*lf
    atoms['K3ml2^2'] = K3m*l2**2; atoms['K3ml2lf'] = K3m*l2*lf; atoms['K3mlf^2'] = K3m*lf**2
    INV = ['z2l2^3','z2l2^2lf','z2l2lf^2','z2lf^3','z3l2^2','z3l2lf','z3lf^2','z2^2l2','z2^2lf',
           'z5','z2z3','z2L3m','K5m','K5q','K4ml2','K4mlf','K4ql2','K4qlf','K3ml2^2','K3ml2lf','K3mlf^2']
    logval = {}
    for nm, (p, q) in LOGPART.items():
        # parse monomial name
        import re
        a = re.search(r'l2\^?(\d*)', nm); b = re.search(r'lf\^?(\d*)', nm)
        ea = int(a.group(1)) if a and a.group(1) else (1 if a else 0)
        eb = int(b.group(1)) if b and b.group(1) else (1 if b else 0)
        logval[nm] = (mp.mpf(p) + mp.mpf(q)*mp.sqrt(5)) * l2**ea * lf**eb
    return atoms, INV, logval

def main():
    D1, D2 = int(sys.argv[1]), int(sys.argv[2])
    here = os.path.dirname(os.path.abspath(__file__))
    W5 = pickle.load(open(os.path.join(here, 'wparts_dec_d520.pkl'),'rb'))
    W8 = pickle.load(open(os.path.join(here, 'wparts_dec_d800.pkl'),'rb'))
    results = []
    for D, W in [(D1, W5), (D2, W8)]:
        atoms, INV, logval = setup(D)
        mp.mp.dps = D + 40
        T = mp.mpf(W['Wt5']) - sum(logval.values())
        at = [(nm, atoms[nm]) for nm in INV]
        cands = fit(T, at, cexp=int(D*0.92), nbest=4)
        print(f"--- R fit at {D} digits (21 invisible atoms) ---")
        report(f"R@{D}", cands)
        results.append(cands)
    r0, r1 = results[0][0], results[1][0]
    same = r0[2] == r1[2]
    print(f">>> R: identical at {D1} and {D2}: {same}")
    if same:
        print("    residuals: %s / %s ; height %s" % (mp.nstr(mp.mpf(r0[1]),3), mp.nstr(mp.mpf(r1[1]),3), r0[0]))
        print(json.dumps({k: (str(v[0]), str(v[1])) for k, v in sorted(r0[2].items())}, indent=1))
        pickle.dump({'rel': r0[2], 'height': r0[0]}, open(os.path.join(here, 'fit_R.pkl'),'wb'))

if __name__ == "__main__":
    main()
'''
open('step_48_rfit.py','w').write(src)
import subprocess, time
PY = '<project>/.venv/bin/python'
t0=time.time()
r = subprocess.run(['nice', PY, 'step_48_rfit.py', '450', '700'], capture_output=True, text=True, cwd='.')
print(r.stdout[-5000:])
if r.stderr: print("STDERR:", r.stderr[-600:])
print("wall %.0fs" % (time.time()-t0))
