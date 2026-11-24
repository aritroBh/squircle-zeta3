import sys, json, math, time
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
def canon_key(parts): return tuple(sorted(parts))
GEN = {'lE13': (Ff(3,2),Ff(1,2)), 'ls13': (Ff(0),Ff(1)), 'l2': (Ff(2),Ff(0)), 'lp3': (Ff(5,2),Ff(1,2)), 'lp3c': (Ff(5,2),Ff(-1,2))}
ca=[Ff(c['x'][0]) for c in J['cands']]; cb=[Ff(c['x'][1]) for c in J['cands']]

def save_json(v, path, extra):
    e=expand(v)
    T_coef = -Ff(e['a_T'])*13**5
    z5_coef = Ff(e['a_z5'])
    li5={f"{ca[i]}+{cb[i]}*sqrt13": Ff(c) for i,c in e['li5'].items()}
    prod={}
    nraw=0
    for k,c in e['prod'].items():
        nraw+=1
        ck='*'.join(canon_key(k))
        prod[ck]=prod.get(ck,Ff(0))+Ff(c)
    prod={k:v for k,v in prod.items() if v}
    out={'D':13, 'T_coef':str(T_coef), 'z5_coef':str(z5_coef),
         'li5':{k:str(v) for k,v in li5.items()},
         'products':{k:str(v) for k,v in prod.items()},
         'letters':{n:f'{a}+{b}*sqrt13' for n,(a,b) in GEN.items()},
         'convention':'T_coef*L(chi13,5)/sqrt13 + z5_coef*zeta(5) + sum li5 + sum products = 0 (same side as d17_final.json)'}
    out.update(extra)
    json.dump(out, open(path,'w'), indent=1)
    coefs=list(li5.values())+list(prod.values())+[z5_coef]
    mn=max(abs(c.numerator) for c in coefs); lm=1
    for c in coefs: lm=math.lcm(lm,c.denominator)
    print(f'{path}: T_coef={T_coef} (={T_coef/13**5}*13^5) z5={z5_coef}; li5 {len(li5)} prod {len(prod)} (raw {nraw}); maxnum {mn}; lcm {lm}')
    return out

vatom=[int(x) for x in json.load(open('out/d13/gcd_vec_atom.json'))['v_atom']]
vnum=[int(x) for x in json.load(open('out/d13/gcd_vec_num.json'))['v_num']]
save_json(vatom, 'd13_gcd.json', {'note':'T-coef = gcd = 36771482407; atom-space max-coef minimised (LLL+Babai+enum); NUMERICAL (LLL), not a proof'})
save_json(vnum, 'd13_gcd_minnum.json', {'note':'T-coef = gcd = 36771482407; expanded max-numerator minimised over T-free (exact-zero) lattice; NUMERICAL (LLL), not a proof'})

# quick resid sanity with correct 13^5 scaling at 300 digits
cands=[dict(x=QD(Ff(c['x'][0]),Ff(c['x'][1]),13), vx=c['vx'], v1m=c['v1m'], h=c['h']) for c in J['cands']]
F=W.make_field(13, split=[3], inert=[2])
for tag,v in (('v_atom',vatom),('v_num',vnum)):
    ee=expand(v); ee['a_T']=Ff(ee['a_T'])*13**5
    r=V.explicit_resid(13, {'D':13,'names':names,'letters':F['letters']}, cands, ee, 300)
    print(tag,'resid@300 (13^5-scaled):', mp.nstr(r,3))
with open('progress.log','a') as f:
    f.write('step_01_gcd: T-coef gcd 36771482407; v_atom atom-max 1550285430850157040 (maxnum 13253867932032456368); v_num maxnum 9070223133044173680 (lcm 1); symbol EXACT 0; factor 143/143 OK; d13_gcd.json + d13_gcd_minnum.json saved\n')
print('progress.log updated')
