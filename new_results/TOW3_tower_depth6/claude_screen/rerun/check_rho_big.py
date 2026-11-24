"""check_rho_big.py -- extend the rho(5), rho(6) check beyond the pickle range (d>52) using TOWN1's truncated engine
(te.rho_level, same engine as the pickle, so this is a range extension, not an independent derivation)."""
import sys, time
sys.path.insert(0, '.')
sys.path.insert(0, '<project>/research/08_theory/TOWN1_tower_large_n')
import trunc_engine as te
import check_data as CD
a5, a6 = CD.load_md_closed_forms()
ev = {5: CD.make_eval(a5, (CD.d_, CD.m_)), 6: CD.make_eval(a6, (CD.d_, CD.m_))}
for d in [int(x) for x in sys.argv[1:]]:
    t0 = time.time()
    S2 = te.stirling_table(d+1)
    R = te.rho_level(d, 6, S2)
    tot = 0; bad = 0
    for (p, q), vals in R.items():
        r0 = CD.Q(vals[0].p0, vals[0].p1, vals[0].p2, vals[0].p3)
        for dl in (5, 6):
            x = vals[dl]
            rr = CD.Q(x.p0, x.p1, x.p2, x.p3) * r0.inv()
            er, ei = ev[dl](d, p)
            tot += 1
            if not (rr.in_Qi() and rr.a == er and rr.c == ei):
                bad += 1; print('BAD', d, p, q, dl)
    print(f'd={d}: {tot} cells (all keys, delta=5,6), bad={bad}  ({round(time.time()-t0,1)}s)', flush=True)
