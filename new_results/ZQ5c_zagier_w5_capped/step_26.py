import os
cwd = os.getcwd()
driver = r'''import sys, json, time, math
from fractions import Fraction as Ff
from collections import defaultdict
import mpmath as mp
sys.path.insert(0,'.')
import w5blib as W
from w5blib import QD
import w3lib_c as X

def fqd(x, D):
    if x.b==0: return f"{x.a}"
    sgn="+" if x.b>0 else "-"
    if x.a==0: return f"{x.b}*sqrt{D}"
    return f"{x.a}{sgn}{abs(x.b)}*sqrt{D}"

def run_D(D, F, cands, tag):
    print(f"===== {tag}: N={len(cands)} =====", flush=True)
    r500 = X.w3_run(D, F, cands, 500)
    r600 = X.w3_run(D, F, cands, 600)
    names = r600['names']; d = r600['d']
    # all relations, expanded, fresh-verified at 700
    forms=[]
    for (P, rr) in ((500,r500),(600,r600)):
        for a in rr['rels']:
            e = X.expand_a(D, F, cands, rr, a)
            res = X.explicit_resid(D, F, cands, e, 700)
            sz = X.symbol_check(cands, e, d)
            nterms = len(e['li3'])+len(e['prod'])+len(e['pi2'])+(1 if e['a_z3'] else 0)
            allc = list(e['li3'].values())+list(e['prod'].values())+list(e['pi2'].values())+[e['a_z3']]
            mxnum = max((abs(x.numerator) for x in allc), default=0)
            forms.append(dict(P=P, e=e, resid=res, symzero=(len(sz)==0), nterms=nterms, mxnum=mxnum, a_T=e['a_T']))
            print(f"  P={P} rel: a_T={e['a_T']} nterms={nterms} maxnum={mxnum} symzero={len(sz)==0} fresh700={mp.nstr(res,3)}", flush=True)
    # agreement: all forms valid (|resid| tiny) and same a_T
    valid = [f for f in forms if abs(f['resid']) < mp.mpf('1e-650') and f['symzero']]
    print(f"{tag}: {len(valid)}/{len(forms)} forms valid @700d; distinct a_T: {sorted(set(f['a_T'] for f in valid))}", flush=True)
    # pick min nterms then min maxnum
    valid.sort(key=lambda f:(f['nterms'], f['mxnum']))
    best = valid[0]
    e = best['e']
    print(f"{tag}: chosen form: P={best['P']} a_T={e['a_T']} nterms={best['nterms']} maxnum={best['mxnum']}", flush=True)
    # pairwise difference exact-zero check (first two valid forms)
    if len(valid)>=2:
        d1, d2 = valid[0]['e'], valid[1]['e']
        dd = dict(a_T=d1['a_T']-d2['a_T'], a_z3=d1['a_z3']-d2['a_z3'])
        li3 = defaultdict(Ff)
        for i,cf in d1['li3'].items(): li3[i]+=cf
        for i,cf in d2['li3'].items(): li3[i]-=cf
        prod = defaultdict(Ff)
        for m,cf in d1['prod'].items(): prod[m]+=cf
        for m,cf in d2['prod'].items(): prod[m]-=cf
        pi2 = defaultdict(Ff)
        for i,cf in d1['pi2'].items(): pi2[i]+=cf
        for i,cf in d2['pi2'].items(): pi2[i]-=cf
        de = dict(a_T=dd['a_T'], a_z3=dd['a_z3'],
                  li3={i:c for i,c in li3.items() if c}, prod={m:c for m,c in prod.items() if c},
                  pi2={i:c for i,c in pi2.items() if c})
        rd = X.explicit_resid(D, F, cands, de, 700)
        print(f"{tag}: difference of two valid forms @700d: {mp.nstr(rd,3)} (exact-zero confirmation)", flush=True)
    allc = list(e['li3'].values())+list(e['prod'].values())+list(e['pi2'].values())+[e['a_z3']]
    lc=1
    for x in allc: lc=lc*x.denominator//math.gcd(lc,x.denominator)
    li5out={fqd(cands[i]['x'], D): str(cf) for i,cf in sorted(e['li3'].items())}
    prodout={}
    for ms,cf in e['prod'].items(): prodout['*'.join(names[ii] for ii in ms)]=str(cf)
    for ii,cf in e['pi2'].items(): prodout[f"pi2*{names[ii]}"]=str(cf)
    out=dict(D=D, weight=3, T_coef=e['a_T'], z3_coef=str(e['a_z3']),
             li3=li5out, products=prodout, names=names,
             letters=[W.fmt_qd(g) for g in F['letters']],
             stats=dict(li3_terms=len(li5out), product_terms=len(prodout),
                        max_numerator=best['mxnum'], lcm_den=lc, total_terms=best['nterms']),
             symbol_zero=True, n_valid_forms=len(valid),
             resid_fresh_700=mp.nstr(best['resid'],4),
             note="NUMERICAL (LLL), not a proof; multiple valid forms differ by exact-zero atom relations")
    json.dump(out, open(f'w3_D{D}.json','w'), indent=1)
    print(f"{tag}: wrote w3_D{D}.json: li3 {len(li5out)}, products {len(prodout)}, maxnum {best['mxnum']}, lcm {lc}", flush=True)
    return out

F5 = W.make_field(5, split=[], inert=[2])
c5 = W.enum_exceptional_fast(F5, B=8)
print(f"D=5: {len(c5)} exceptional, letters {F5['names']}", flush=True)
run_D(5, F5, c5, 'D5w3')
with open('progress.log','a') as f: f.write("step_04_w3_d5: w3_D5.json (verified multi-precision)\n")

F13 = W.make_field(13, split=[3], inert=[2])
c13 = W.enum_exceptional_fast(F13, B=6)
c13 = [c for c in c13 if c['h']<=3]
print(f"D=13: {len(c13)} exceptional (h<=3, B=6)", flush=True)
run_D(13, F13, c13, 'D13w3')
with open('progress.log','a') as f: f.write("step_04_w3_d13: w3_D13.json (verified multi-precision)\n")
print("W3 V2 DONE", flush=True)
'''
open('step_04_w3v2.py','w').write(driver)
import subprocess, time
t0=time.time()
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_04_w3v2.py'],
                   capture_output=True, text=True, cwd=cwd, timeout=3000)
print("elapsed", time.time()-t0)
print(r.stdout)
print("STDERR:", r.stderr[-2000:] if r.stderr else "")
