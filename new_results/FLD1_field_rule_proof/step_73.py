WD='<project>/research/08_theory/FLD1_field_rule_proof'
with open(WD+'/progress.log','a') as fh:
    fh.write("step31: k=6 symbolic mu == table (n<=60); to_zsum; PROPER (deg-2), weight budget (<=7-wt(z), saturated), poles -(5/2)..1/2. 'PROVED FOR k<=6' for structural items\n")
code = r'''"""step32_k6construct.py -- full Lemma C construction + all-n proof for k=6 (the algorithm scales)."""
import sys, time, pickle
from fractions import Fraction as Fr
sys.path.insert(0,'.')
src=open("step_25_construct.py").read(); ns={}
exec(src[:src.index("if __name__")], ns)
construct=ns['construct']; symmetrise=ns['symmetrise']; correct=ns['correct']
even_form=ns['even_form']; odd_form=ns['odd_form']; toK=ns['toK']
import zsum
from zsum import lin
Zin = pickle.load(open('mu_zsum_K6.pkl','rb'))  # k=6 only
t0=time.time()
store={}; allok=True
for zkey, E in sorted(Zin.items()):
    Elist=[(w,c) for w,c in E]
    terms=construct(Elist); Fsym=symmetrise(terms); Fsym,bad=correct(Fsym,6)
    store[zkey]=Fsym
    ev=even_form(Fsym); od=odd_form(Fsym)
    muE={tuple(w): toK(c) for w,c in Elist}
    diff=lin((1,ev),(-1,muE))
    ok = len(diff)==0 and len(od)==0
    allok = allok and ok
    print(f'k=6 z={zkey}: {len(terms)} terms -> {len(Fsym)} (corr {bad}); even==mu: {len(diff)==0}; odd==0: {len(od)==0}  [{time.time()-t0:.0f}s]', flush=True)
pickle.dump(store, open('hpl_constructed_K6.pkl','wb'))
print('k=6 ALL-N construction proof (even==mu & odd==0, all zeta-monomials):', allok)
'''
open(WD+'/step_32_k6construct.py','w').write(code)
import subprocess, time
VENV='<project>/.venv/bin/python'
t=time.time()
r = subprocess.run(['nice',VENV,'step_32_k6construct.py'],capture_output=True,text=True,cwd=WD, timeout=1700)
print("wall",round(time.time()-t,1),"s")
for line in r.stdout.splitlines():
    if 'residual' in line: continue
    print(line)
print(r.stderr[-2000:] if r.returncode else '')