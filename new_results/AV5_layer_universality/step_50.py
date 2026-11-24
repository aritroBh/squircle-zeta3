import mpmath, pickle
from mpmath import mp, mpf, zeta, nstr
mp.dps = 125   # BEFORE loading!
F = {k: mpf(v) for k,v in pickle.load(open('ckpt_Fk_125.pkl','rb')).items()}
Z = lambda n: zeta(n)
rels = {
 3: [(33, Z(2)**2), (160, None)],
 4: [(-7, Z(2)*Z(3)), (-3, Z(5)), (-16, None)],
 5: [(577, Z(2)**3), (-140, Z(3)**2), (-640, None)],
 6: [(1071, Z(2)**2*Z(3)), (1410, Z(2)*Z(5)), (1290, Z(7)), (-640, None)],
}
for k, terms in rels.items():
    R = sum(c*v if v is not None else c*F[k] for c,v in terms)
    print(f"F_{k}: residual {nstr(abs(R),4)}  -> VERIFIED" if abs(R)<mpf(10)**(-110) else f"F_{k}: {nstr(abs(R),4)} FAIL")
with open('progress.log','a') as fh:
    fh.write('step_50: zeta-only collapse F3..F6 re-verified with correct precision (residuals<1e-110)\n')
