"""step_05_weight2_values.py: 64 weight-2 word values + quotient basis at ~325 digits.
NOTE: check_conversion bakes roots/N at import precision; recompute after raising ctx.prec."""
import sys, pickle, time, json
sys.dont_write_bytecode = True
sys.path.insert(0, 'Q')
from flint import ctx, acb, arb
import check_conversion as CC
ctx.prec = 1100
ctx.cap = 2048
CC.N = 2048
CC.roots = [acb(0, arb.pi()*j/4).exp() for j in range(8)]
CC.roots[0], CC.roots[2], CC.roots[4], CC.roots[6] = acb(1), acb(0,1), acb(-1), acb(0,-1)
from check_conversion import h
import mpmath as mp
mp.mp.dps = 335
from relations import words

def s330(z):
    return mp.mpf(z.mid().str(340, radius=False, more=True))

def main():
    cols = words(2)
    t0 = time.monotonic()
    vals = {w: h(w) for w in cols}
    print(f'64 values computed ({time.monotonic()-t0:.0f}s)', flush=True)
    worst = max(max(v.real.rad(), v.imag.rad()) for v in vals.values())
    print('max ball radius:', mp.nstr(mp.mpf(worst.str(30, radius=False, more=True)), 8), flush=True)
    with open('audit_weight2.pickle','rb') as f:
        data = pickle.load(f)
    from flint import fmpq_mat
    rows = data['rows']
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
    out = {'prec': 1100, 'cap': 1150, 'cols': [list(w) for w in cols],
           'free_words': [list(cols[j]) for j in free],
           'values': {str(list(w)): [mp.nstr(s330(vals[w].real), 330),
                                     mp.nstr(s330(vals[w].imag), 330)] for w in cols}}
    with open('weight2_values.json','w') as f:
        json.dump(out, f)
    for j in free:
        w = cols[j]
        print(w, mp.nstr(mp.mpc(s330(vals[w].real), s330(vals[w].imag)), 25), flush=True)
    # spot check a known value: h((4,4)) should be related to Li2(-1) = -pi^2/12
    import mpmath as mm
    print('check h(4,4) vs Li2(-1)=-pi^2/12:',
          mp.nstr(mp.mpc(s330(vals[(4,4)].real), s330(vals[(4,4)].imag)) + mp.pi**2/12, 10))

main()
