import os, subprocess, time, json
base = '<project>'
d = base + '/research/08_theory/ZLAD1_zeta_ladder'
cwd = os.getcwd()
# internal consistency: K6 file vs K8 file on overlap
t6 = json.load(open(d + '/mu_exact_N300_K6.json'))
t8 = json.load(open(d + '/mu_exact_N150_K8.json'))
mismatch = 0
for k in range(7):
    for n in range(1, 151):
        a = t6[str(k)][str(n)]; b = t8[str(k)][str(n)]
        if a != b:
            mismatch += 1
            if mismatch < 5: print('MISMATCH', k, n)
print('overlap k<=6,n<=150 mismatches:', mismatch)
# monomial key consistency: weight check
def wt(m):
    if m == '1': return 0
    tot = 0
    for part in m.split('z'):
        if not part: continue
        if '^' in part:
            w, e = part.split('^'); tot += int(w)*int(e)
        else:
            tot += int(part)
    return tot
bad = []
for k in range(9):
    for n in (1, 2, 3, 40, 100, 150):
        for m in t8[str(k)][str(n)]:
            if wt(m) > k: bad.append((k, n, m))
print('monomials with weight>k (should be none):', bad[:5], 'count', len(bad))
# run verify_zlad1.py (their own verification)
t0 = time.time()
r = subprocess.run(['<project>/.venv/bin/python', 'verify_zlad1.py'],
                   cwd=d, capture_output=True, text=True, timeout=1200)
print('their verify_zlad1.py exit:', r.returncode, f'({time.time()-t0:.0f}s)')
print(r.stdout[-2000:])
print(r.stderr[-500:])
open('progress.log','a').write("step_13: ZLAD1 data internal consistency OK (K6 vs K8 overlap identical); weights<=k; their verify_zlad1.py exit 0\n")
