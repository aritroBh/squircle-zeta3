import sys, time, pickle
sys.path.insert(0,'.')
from oct9c_common import *
D = load_all()
Gargs, Dargs = rhs_terms_235()
sols7, negw7, logprods7 = rhs_terms_2357()

def remainder(e, dps):
    mp.mp.dps = dps
    eng = get_engine(e)
    NUM = atom_values(eng, e, D['lam'], dps)
    if e == 0:
        W3 = D['Wp0']['W3_e0']; sol = D['s0']['solB0']
        RHS = eval_rhs_2357(sol, sols7, negw7, logprods7, dps)
    else:
        W3 = D['s123'][e]['W3']; sol = D['s123'][e]['sol']
        RHS = eval_rhs_235(sol, Gargs, Dargs, dps)
    W3val = mp.mpf(sp.N(W3.subs(list(NUM.items())), dps))
    return W3val - RHS, W3val, RHS

out = {}
for dps in (250, 340):
    mp.mp.dps = dps
    z2 = mp.pi**2/6
    for e in range(4):
        t0=time.time()
        R, W3val, RHS = remainder(e, dps)
        inv = [mp.zeta(3), z2*mp.log(2), z2*mp.log(3), z2*mp.log(5)] + ([z2*mp.log(7)] if e==0 else [])
        rels = mp.pslq([R]+inv, maxcoeff=10**15, maxsteps=10000)
        err = None
        if rels:
            err = abs(rels[0]*R + mp.fsum([mp.mpf(r_)*v for r_,v in zip(rels[1:], inv)]))
        out[(e,dps)] = {'R': mp.nstr(R, dps-5), 'rels': rels, 'err': mp.nstr(err,3) if err is not None else None}
        print(f"edge {e} dps={dps}: R={mp.nstr(R,12)} rels={rels} err={mp.nstr(err,3) if err else None} [{time.time()-t0:.0f}s]", flush=True)
pickle.dump(out, open('task1_pslq_raw.pkl','wb'))
print("saved task1_pslq_raw.pkl")
