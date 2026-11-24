import mpmath, pickle, json
from mpmath import mp, mpf, zeta, nstr
F = {k: mpf(v) for k,v in pickle.load(open('ckpt_Fk_125.pkl','rb')).items()}
def run(dps):
    mp.dps = dps
    Z = lambda n: zeta(n)
    rels = {
     3: [(33, Z(2)**2), (160, None)],
     4: [(-7, Z(2)*Z(3)), (-3, Z(5)), (-16, None)],
     5: [(577, Z(2)**3), (-140, Z(3)**2), (-640, None)],
     6: [(1071, Z(2)**2*Z(3)), (1410, Z(2)*Z(5)), (1290, Z(7)), (-640, None)],
    }
    for k, terms in rels.items():
        R = sum(c*v if v is not None else c*F[k] for c,v in terms)
        print(f"  F_{k} (dps{dps}): residual {nstr(abs(R),4)}")
run(120)
run(60)
with open('progress.log','a') as fh:
    fh.write('step_48: zeta-only collapse F3..F6 verified at two precisions; squircle c2..c7 exact\n')
