import sys, json, time
from fractions import Fraction as Ff
import mpmath as mp
sys.path.insert(0,'.')
import w5blib as W, w5verify as V
from w5blib import QD
gl = json.load(open('out/d13/genuine_lattice_cap.json'))
atom_names = gl['atom_names']
J = json.load(open('out/d13/candidates_cap.json'))
Kk = json.load(open('out/d13/kernel5_cap.json'))
K=Kk['K']; comps=Kk['comps']; names=J['names']
keepc=[0,1,2,3]+list(range(23,51))
def expand(v): return V.expand_a(13, J['cands'], K, comps, keepc, atom_names, v, names)
cands=[dict(x=QD(Ff(c['x'][0]),Ff(c['x'][1]),13), vx=c['vx'], v1m=c['v1m'], h=c['h']) for c in J['cands']]
F=W.make_field(13, split=[3], inert=[2])
vatom=[int(x) for x in json.load(open('out/d13/gcd_vec_atom.json'))['v_atom']]
vnum=[int(x) for x in json.load(open('out/d13/gcd_vec_num.json'))['v_num']]
for tag,v in (('v_atom',vatom),('v_num',vnum)):
    ee=expand(v)
    t0=time.time(); sym=V.symbol_check5(cands, ee, names)
    t1=time.time(); nchk,bad=V.factorization_check(cands, ee, F)
    t2=time.time(); r300=V.explicit_resid(13, {'D':13,'names':names,'letters':F['letters']}, cands, ee, 300)
    print(f'{tag}: symbol nonzero={len(sym)} ({t1-t0:.0f}s); factor checked {nchk} bad {bad} ({t2-t1:.0f}s); resid@300 {mp.nstr(r300,3)} ({time.time()-t2:.0f}s)')
