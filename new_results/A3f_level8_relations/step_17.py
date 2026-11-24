open('progress.log','a').write(
 'step_02: AUDIT PASS w1-3 exact Q quotients 3/10/37 (fmpq rref + 3 primes).\n'
 'step_03: AUDIT PASS all w1-3 relation rows numerically zero, worst 4.8e-58 (acb prec=200).\n'
 'step_04: AUDIT PASS w4 quotient 142 mod second prime 1000000021; all 17625 rows numerically zero (worst 2.9e-57); D residue 141/142.\n')

code = r'''"""step_05_weight2_values.py: 64 weight-2 word values + quotient basis at ~330 digits."""
import sys, pickle, time, json
sys.dont_write_bytecode = True
sys.path.insert(0, 'Q')
from flint import ctx
import check_conversion as CC
# raise precision BEFORE any cached evaluation
ctx.prec = 1100          # ~331 digits
ctx.cap = 1150           # series truncation ~ 2^-1150 ~ 1e-346
from check_conversion import h, approx
import mpmath as mp
mp.mp.dps = 330
from relations import words

def main():
    cols = words(2)
    t0 = time.monotonic()
    vals = {}
    for w in cols:
        v = h(w)
        vals[w] = v
    print(f'64 values computed ({time.monotonic()-t0:.0f}s)', flush=True)
    # quotient basis from audit pickle (exact Q rref pivot structure)
    with open('audit_weight2.pickle','rb') as f:
        data = pickle.load(f)
    from flint import fmpq_mat
    rows, cols2 = data['rows'], data['cols']
    assert cols2 == cols
    ids = {w:i for i,w in enumerate(cols)}
    M = fmpq_mat(len(rows), len(cols))
    for i,row in enumerate(rows):
        for w,c in row: M[i, ids[w]] = int(c)
    rr, rank = M.rref(inplace=True)
    pp = []; start=0
    for i in range(rank):
        while not rr[i,start]: start+=1
        pp.append(start); start+=1
    free = [j for j in range(len(cols)) if j not in set(pp)]
    print('rank', rank, 'free words:', [cols[j] for j in free], flush=True)
    out = {'prec': 1100, 'cap': 1150, 'free_words': [list(cols[j]) for j in free],
           'values': {str(list(w)): [mp.nstr(mp.mpf(vals[w].real.mid().str(340, radius=False, more=True)), 335),
                                     mp.nstr(mp.mpf(vals[w].imag.mid().str(340, radius=False, more=True)), 335)] for w in cols}}
    with open('weight2_values.json','w') as f:
        json.dump(out, f)
    for j in free:
        w = cols[j]
        print(w, mp.nstr(mp.mpc(mp.mpf(vals[w].real.mid().str(340,radius=False,more=True)),
                                mp.mpf(vals[w].imag.mid().str(340,radius=False,more=True))), 30), flush=True)
    # sanity: radius check
    worst = max(max(vals[w].real.rad(), vals[w].imag.rad()) for w in cols)
    print('max ball radius:', mp.nstr(mp.mpf(worst.str(30, radius=False, more=True)), 8))

main()
'''
open('step_05_weight2_values.py','w').write(code)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python','-B','step_05_weight2_values.py'],
                   capture_output=True, text=True, timeout=2400)
print(r.stdout); print('STDERR:', r.stderr[-2000:])
