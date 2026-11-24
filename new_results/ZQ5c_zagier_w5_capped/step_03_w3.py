import sys, json, time, math
from fractions import Fraction as Ff
from collections import defaultdict
import mpmath as mp
sys.path.insert(0,'.')
import w5blib as W
from w5blib import QD
from w3lib_c import w3_rel, expand

def run_D(D, F, cands, tag):
    print(f"===== {tag}: N={len(cands)} =====", flush=True)
    r500 = w3_rel(D, F, cands, 500)
    r600 = w3_rel(D, F, cands, 600)
    if r500['best'] is None or r600['best'] is None:
        print(f"{tag}: NO RELATION at one precision"); return None
    e500 = expand(D, F, cands, r500)
    e600 = expand(D, F, cands, r600)
    # compare relations (normalize by T_coef)
    def norm(e):
        t=e['a_T']
        return (t, e['a_z3'],
                {k: v/t for k,v in e['li3'].items()},
                {k: v/t for k,v in e['prod'].items()},
                {k: v/t for k,v in e['pi2'].items()})
    n5, n6 = norm(e500), norm(e600)
    same = (n5==n6)
    print(f"{tag}: two-precision agreement: {same}; T_coef {e500['a_T']} (P500), {e600['a_T']} (P600)", flush=True)
    if not same:
        print("  MISMATCH -- keeping 600 version")
    e = e600; rr = r600
    # exact symbol check (weight 3): sum li3c*sym_li3 + sum prodc*prod_symbol = 0
    d = rr['d']; names = rr['names']
    T=defaultdict(Ff)
    for i,cf in e['li3'].items():
        s=W.symbol_liw_dense(cands[i]['vx'],cands[i]['v1m'],d,3)
        for idx,v in enumerate(s):
            if v: T[idx]+=cf*v
    for ms,cf in e['prod'].items():
        v=W.prod_symbol_dense(list(ms),d,3)
        for i2,vv in enumerate(v):
            if vv: T[i2]+=cf*vv
    Tz={k:v for k,v in T.items() if v}
    symzero = (len(Tz)==0)
    print(f"{tag}: exact symbol of expanded form: {'ZERO' if symzero else f'NONZERO {len(Tz)}'}", flush=True)
    # fresh residual at 700 digits with the EXPLICIT terms
    P=700; mp.mp.dps=P+40
    sD=mp.sqrt(D)
    lg=[mp.log(abs(g.emb(P+40))) for g in F['letters']]
    L3=W.Lchi_hurwitz(D,3,P)
    tot = Ff(e['a_T'])*(L3/sD) - Ff(e['a_z3'])*mp.zeta(3)
    for i,cf in e['li3'].items():
        tot -= cf*W.li3(cands[i]['x'].emb(P+40), P+20)
    for ms,cf in e['prod'].items():
        pv=mp.mpf(1)
        for ii in ms: pv*=lg[ii]
        tot -= cf*pv
    for ii,cf in e['pi2'].items():
        tot -= cf*mp.pi**2*lg[ii]
    resid700 = tot
    print(f"{tag}: fresh residual of expanded form @700d: {mp.nstr(resid700,4)}", flush=True)
    # stats
    allc = list(e['li3'].values())+list(e['prod'].values())+list(e['pi2'].values())+[Ff(e['a_z3'])]
    mxnum = max(abs(x.numerator) for x in allc) if allc else 0
    lc=1
    for x in allc:
        lc=lc*x.denominator//math.gcd(lc,x.denominator)
    # JSON
    def fqd(x):
        if x.b==0: return f"{x.a}"
        sgn="+" if x.b>0 else "-"
        return f"{x.a}{sgn}{abs(x.b)}*sqrt{D}" if x.a!=0 else (f"{x.b}*sqrt{D}")
    li5out={fqd(cands[i]['x']): str(cf) for i,cf in sorted(e['li3'].items())}
    prodout={}
    for ms,cf in e['prod'].items():
        prodout['*'.join(names[ii] for ii in ms)]=str(cf)
    for ii,cf in e['pi2'].items():
        prodout[f"pi2*{names[ii]}"]=str(cf)
    out=dict(D=D, weight=3, T_coef=e['a_T'], z3_coef=str(e['a_z3']),
             li3=li5out, products=prodout,
             letters=[f"log|{g}|=log({fmt})" for g,fmt in zip(F['letters'],[W.fmt_qd(g) for g in F['letters']])],
             names=names,
             stats=dict(li3_terms=len(li5out), product_terms=len(prodout),
                        max_numerator=mxnum, lcm_den=lc),
             symbol_zero=bool(symzero),
             agreement_500_600=bool(same),
             resid_P500=mp.nstr(e500['resid'],4), resid_P600=mp.nstr(e600['resid'],4),
             resid_fresh_700=mp.nstr(resid700,4),
             note="NUMERICAL (LLL), not a proof")
    json.dump(out, open(f'w3_D{D}.json','w'), indent=1)
    print(f"{tag}: wrote w3_D{D}.json: li3 {len(li5out)} terms, products {len(prodout)}, maxnum {mxnum}, lcm {lc}", flush=True)
    return out

# D=5
F5 = W.make_field(5, split=[], inert=[2])
c5 = W.enum_exceptional_fast(F5, B=8)
print(f"D=5: {len(c5)} exceptional, letters {F5['names']}", flush=True)
run_D(5, F5, c5, 'D5w3')
with open('progress.log','a') as f: f.write("step_03_w3_d5: w3_D5.json done\n")

# D=13
F13 = W.make_field(13, split=[3], inert=[2])
c13 = W.enum_exceptional_fast(F13, B=6)
c13 = [c for c in c13 if c['h']<=3]
print(f"D=13: {len(c13)} exceptional (h<=3, B=6), letters {F13['names']}", flush=True)
run_D(13, F13, c13, 'D13w3')
with open('progress.log','a') as f: f.write("step_03_w3_d13: w3_D13.json done\n")
print("W3 ALL DONE", flush=True)
