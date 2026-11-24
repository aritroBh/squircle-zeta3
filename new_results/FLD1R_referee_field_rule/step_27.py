import subprocess, json, sys
r = subprocess.run(['nice','<project>/.venv/bin/python','step_21_fastmu.py','7','30'],
                   capture_output=True, text=True, cwd='.')
print(r.stdout[-700:], r.stderr[-300:] if r.stderr else '')
# compare the two tables
mine = json.load(open('mu7_indep_N30.json'))
ref  = json.load(open('mu_exact_K7_N30.json'))
ok=True
for k in range(8):
    for nv in range(1,31):
        A = mine[str(k)][str(nv)]
        B = {kk.strip('()'):vv for kk,vv in ref[str(k)][str(nv)].items()}
        keys = set(A)|set(B)
        for kk in keys:
            from fractions import Fraction as Fr
            if Fr(A.get(kk,'0')) != Fr(B.get(kk,'0')):
                ok=False; print('MISMATCH', k, nv, kk, A.get(kk), B.get(kk))
print('INDEPENDENT (sympy series) == FLD1 step_21 (Stirling) for k<=7, n<=30:', ok)
open('progress.log','a').write(f'step03b: indep mu_{{k,n}} (k<=7,n<=30) == FLD1 step_21 Stirling table: {ok}\n')
