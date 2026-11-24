import sys, os, time
workdir = '<project>/research/08_theory/TOW2R_review'
sys.path.insert(0, workdir)
import tow2r_lib as RL
sys.path.insert(0, os.path.join(workdir, 'tow2_copy'))
import towlib
towlib.WTS = '<project>/research/08_theory/WTS_weight_k_structure/theory'

t0 = time.time()
bad = []; tested = 0
for dd in range(4, 13):
    for key in towlib.keys_of_degree(dd):
        p, q = key
        mine = RL.rho_first_principles(p, q, 4)
        for delta in range(0, 5):
            ref = towlib.rho(key, delta)
            tested += 1
            if mine[delta] != ref: bad.append((key, delta))
print(f'rho_first_principles vs binteg, d=4..12, delta<=4: {tested} cells, {len(bad)} bad  ({time.time()-t0:.1f}s)')
if bad: print(bad[:8])
# calibrate S_poly cost
t0=time.time(); S30 = RL.GI.S_poly(30, 15); print('S_poly(30,15):', round(time.time()-t0,2),'s')
t0=time.time(); S40 = RL.GI.S_poly(40, 20); print('S_poly(40,20):', round(time.time()-t0,2),'s')
