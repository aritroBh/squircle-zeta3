import sys, pickle, time
sys.path.insert(0,'.')
Z6 = pickle.load(open('mu_zsum_K6.pkl','rb'))
print('k=6 words per zeta-mono:', {zk:len(v) for zk,v in Z6.items()})
Z = pickle.load(open('mu_zsum_K7.pkl','rb'))
# reuse step_25's machinery without __main__
src=open('step_25_construct.py').read(); ns={}
exec(src[:src.index("if __name__")], ns)
construct=ns['construct']; symmetrise=ns['symmetrise']; correct=ns['correct']
even_form=ns['even_form']; odd_form=ns['odd_form']; toK=ns['toK']
import zsum
from zsum import lin
def wt(zkey): return sum((i+2)*e for i,e in enumerate(zkey))
store={}
t0=time.time()
todo = [(zk,E) for zk,E in sorted(Z.items()) if zk != (0,0,0,0,0,0)]
for zkey, E in todo:
    tt=time.time()
    terms = construct(E)
    Fsym = symmetrise(terms)
    Fsym, bad = correct(Fsym, 7)
    ev = even_form(Fsym); od = odd_form(Fsym)
    muE = {tuple(w): toK(c) for w,c in E}
    diff = lin((1,ev),(-1,muE))
    store[zkey]=Fsym
    print(f'z={zkey} (wt {wt(zkey)}): {len(terms)} raw -> {len(Fsym)} terms (corr {bad}); even==mu: {len(diff)==0}; odd==0: {len(od)==0}; {time.time()-tt:.0f}s', flush=True)
    if diff: print('   EVEN residual:', list(diff.items())[:2])
    if od: print('   ODD residual:', list(od.items())[:2])
pickle.dump(store, open('hpl_constructed_K7_nonrat.pkl','wb'))
print('non-rational monomials done in', round(time.time()-t0,1),'s; saved hpl_constructed_K7_nonrat.pkl')
open('progress.log','a').write(f'step05a: k=7 construction+all-n identity for 14 non-rational zeta-monomials: done in {time.time()-t0:.0f}s\n')
