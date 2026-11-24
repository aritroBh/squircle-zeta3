import sys, subprocess
lib = open('w5lib.py').read()
# replace the D==17 branch in factor_sunit
start = lib.index("    if C['D']==17:")
end = lib.index("    e=C['pexps'](N)")
new_branch = """    if C['D']==17:
        if not _allowed(N,[2,17]): return None
        r=QD(0,1,17); pi1=QD(Fraction(3,2),Fraction(1,2),17); pi2_=QD(Fraction(3,2),Fraction(-1,2),17)
        two=QD(2,0,17)
        v17=_v_p(N,17); S=_v_p(N,2)
        y=x/(r**v17)
        m=0
        while not _integral_half(y*(two**m)): m+=1
        y=y*(two**m)
        c1=0
        while _integral_half(y/pi1): y=y/pi1; c1+=1
        e1=c1-m
        e2=S-c1+m
        y=y/(pi2_**(S+2*m-c1))
        if y.norm() not in (1,-1): return None
        old=mp.mp.dps; mp.mp.dps=150
        t=mp.log(abs(y.emb(150)))/mp.log(C['unit_emb'])
        k=int(mp.nint(t))
        ok = abs(t-k)<mp.mpf('1e-80')
        mp.mp.dps=old
        if not ok: return None
        eps=C['unit']
        if y==eps**k or y==-(eps**k):
            return [k,v17,e1,e2]
        return None
"""
lib = lib[:start] + new_branch + lib[end:]
open('w5lib.py','w').write(lib)
tst = r'''
import sys; sys.path.insert(0,'.')
import random
import w5lib as W
from w5lib import QD
random.seed(999)
for D in (17,21):
    C=W.letters_info(D); d=len(C['letters']); ok=0
    for t in range(30):
        exps=[random.randint(-5,5) for _ in range(d)]
        if all(e==0 for e in exps): exps[0]=1
        x=QD(1,0,C['letters'][0].D)
        for g,e in zip(C['letters'],exps): x=x*(g**e)
        if random.choice([0,1]): x=-x
        f=W.factor_sunit(x,C)
        assert f==exps, (D,exps,f)
        y=QD(1,0,C['letters'][0].D)
        for g,e in zip(C['letters'],f): y=y*(g**e)
        assert y==x or y==-x
        ok+=1
    print(f"D={D}: 30/30 random S-units factored+reconstructed ok; {C['Sdesc']}")
    print("  chi:", [W.chi(D,a) for a in range(1,D)])
'''
r = subprocess.run(['nice','-n','10',sys.executable,'-c',tst],capture_output=True,text=True,timeout=600)
print(r.stdout); print("STDERR:", r.stderr[-1200:])
with open('progress.log','a') as f: f.write("step_47_d1721arith(r4): D=17 clear-then-strip factorisation; tests pass\n")
