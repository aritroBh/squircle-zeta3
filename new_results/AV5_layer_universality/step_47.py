import mpmath, pickle, json
from mpmath import mp, mpf, zeta, log, sqrt, gamma, nstr
mp.dps = 125
F = {k: mpf(v) for k,v in pickle.load(open('ckpt_Fk_125.pkl','rb')).items()}
Z = lambda n: zeta(n)
cands = {
 3: -sp_Z2 := None,
}
# verify the zeta-only relations by direct residual
rels = {
 3: [(33, Z(2)**2), (160, None)],
 4: [(-7, Z(2)*Z(3)), (-3, Z(5)), (-16, None)],
 5: [(577, Z(2)**3), (-140, Z(3)**2), (-640, None)],
 6: [(1071, Z(2)**2*Z(3)), (1410, Z(2)*Z(5)), (1290, Z(7)), (-640, None)],
}
for k, terms in rels.items():
    R = sum(c*v if v is not None else c*F[k] for c,v in terms)
    print(f"F_{k}: residual {nstr(abs(R),4)}  (data precision ~1e-115)")
# also at dps 60 for two-precision evidence
mp.dps = 60
F60 = {k: mpf(nstr(v,55)) for k,v in F.items()}
for k, terms in rels.items():
    R = sum(c*v if v is not None else c*F60[k] for c,v in terms)
    print(f"F_{k} (dps60): residual {nstr(abs(R),4)}")
