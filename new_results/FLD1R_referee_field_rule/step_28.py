import json
from fractions import Fraction as Fr
# load step_21's machinery without __main__
src = open('step_21_fastmu.py').read()
ns = {}
exec(src[:src.index("if __name__")], ns)
D = ns['moments'](7)
ref = {}
for k in range(8):
    ref[k] = {}
    for nv in range(1, 31):
        P = ns['mu'](k, nv, 7, D)
        ref[k][nv] = {','.join(map(str,e)): str(v) for e,v in P.items()}
mine = json.load(open('mu7_indep_N30.json'))
ok=True
for k in range(8):
    for nv in range(1,31):
        A, B = mine[str(k)][str(nv)], ref[k][nv]
        for kk in set(A)|set(B):
            if Fr(A.get(kk,'0')) != Fr(B.get(kk,'0')):
                ok=False; print('MISMATCH',k,nv,kk,A.get(kk),B.get(kk))
print('INDEPENDENT sympy-series table == FLD1 step_21 Stirling table (k<=7, n<=30):', ok)
json.dump({str(k):{str(nv):ref[k][nv] for nv in ref[k]} for k in ref}, open('mu_exact_K7_N30.json','w'))
print('saved mu_exact_K7_N30.json')
open('progress.log','a').write(f'step03c: indep mu_{{k,n}} (k<=7,n<=30) == FLD1 Stirling table (all 15 zeta-monomials, exact): {ok}\n')
