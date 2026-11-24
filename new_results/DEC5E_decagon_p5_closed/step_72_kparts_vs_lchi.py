"""step_72_kparts_vs_lchi.py - task (3): do the 16 golden Li32/Li41 K-parts reduce
to the 22-atom classical+L(chi5) basis? Two precisions, identical vectors required.
Usage: python step_72_kparts_vs_lchi.py D1 D2
"""
import sys, os, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from q5_lll import fit
HERE = os.path.dirname(os.path.abspath(__file__))

def atoms(D):
    mp.mp.dps = D + 40
    phi = (1+mp.sqrt(5))/2
    l2 = mp.log(2); lf = mp.log(phi); z2 = mp.zeta(2); z3 = mp.zeta(3); z5 = mp.zeta(5)
    L3m = mp.polylog(3, -1/phi); L4m = mp.polylog(4, -1/phi); L4q = mp.polylog(4, 1/phi**2)
    L5m = mp.polylog(5, -1/phi); L5q = mp.polylog(5, 1/phi**2)
    K5m = L5m + lf**5/120; K5q = L5q - mp.mpf(2)/15*lf**5
    K4m = L4m - lf**4/24;  K4q = L4q + lf**4/3
    K3m = L3m + lf**3/6
    A = {}
    for a in range(4):
        b = 3-a
        nm = 'z2' + (('l2^%d'%a if a>1 else 'l2'*a) + ('lf^%d'%b if b>1 else 'lf'*b))
        A[nm] = z2*l2**a*lf**b
    for a in range(3):
        b = 2-a
        nm = 'z3' + (('l2^%d'%a if a>1 else 'l2'*a) + ('lf^%d'%b if b>1 else 'lf'*b))
        A[nm] = z3*l2**a*lf**b
    A['z2^2l2'] = z2**2*l2; A['z2^2lf'] = z2**2*lf
    A['z5'] = z5; A['z2z3'] = z2*z3
    A['z2L3m'] = z2*L3m
    A['K5m'] = K5m; A['K5q'] = K5q
    A['K4ml2'] = K4m*l2; A['K4mlf'] = K4m*lf
    A['K4ql2'] = K4q*l2; A['K4qlf'] = K4q*lf
    A['K3ml2^2'] = K3m*l2**2; A['K3ml2lf'] = K3m*l2*lf; A['K3mlf^2'] = K3m*lf**2
    Lchi = pickle.load(open(os.path.join(HERE, 'Lchi5_1550.pkl'),'rb'))
    A['Lchi5'] = mp.mpf(Lchi['Lchi5_5'])
    INV = ['z2l2^3','z2l2^2lf','z2l2lf^2','z2lf^3','z3l2^2','z3l2lf','z3lf^2','z2^2l2','z2^2lf',
           'z5','z2z3','z2L3m','K5m','K5q','K4ml2','K4mlf','K4ql2','K4qlf','K3ml2^2','K3ml2lf','K3mlf^2','Lchi5']
    return A, INV

def main():
    D1, D2 = int(sys.argv[1]), int(sys.argv[2])
    kpool = pickle.load(open(os.path.join(HERE, 'kpool_1550.pkl'),'rb'))
    n_pass = 0
    for nm in sorted(kpool):
        rels = []
        for D in (D1, D2):
            A, INV = atoms(D)
            mp.mp.dps = D + 40
            T = mp.mpf(kpool[nm]['K'])
            at = [(x, A[x]) for x in INV]
            c = fit(T, at, cexp=int(D*0.92), nbest=1)
            rels.append(c[0])
        same = rels[0][2] == rels[1][2]
        ok = same and mp.log10(mp.mpf(rels[0][1])+mp.mpf(10)**(-D1-60)) < -0.8*D1 and mp.log10(mp.mpf(rels[1][1])+mp.mpf(10)**(-D2-60)) < -0.8*D2
        if ok: n_pass += 1
        print("%-7s identical: %-5s height %.2e  resid %s / %s  %s" % (
            nm, same, rels[0][0], mp.nstr(mp.mpf(rels[0][1]),2), mp.nstr(mp.mpf(rels[1][1]),2),
            "REDUCES" if ok else "*** DOES NOT REDUCE ***"))
    print("PASS: %d / 16 K-parts reduce to the 22-atom basis" % n_pass)

if __name__ == "__main__":
    main()
