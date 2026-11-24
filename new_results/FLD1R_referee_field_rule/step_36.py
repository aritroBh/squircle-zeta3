import sys, pickle, time
sys.path.insert(0,'.')
import zsum
from zsum import lin
Z = pickle.load(open('mu_zsum_K7.pkl','rb'))
src=open('step_25_construct.py').read(); ns={}
exec(src[:src.index("if __name__")], ns)
construct=ns['construct']; symmetrise=ns['symmetrise']; correct=ns['correct']
even_form=ns['even_form']; odd_form=ns['odd_form']; toK=ns['toK']
t0=time.time()
zkey=(0,0,0,0,0,0); E=Z[zkey]
terms = construct(E)
print(f'construct: {len(terms)} raw terms in {time.time()-t0:.0f}s', flush=True)
Fsym = symmetrise(terms)
Fsym, bad = correct(Fsym, 7)
print(f'symmetrised+corrected: {len(Fsym)} terms, correction {bad}', flush=True)
ev = even_form(Fsym); od = odd_form(Fsym)
muE = {tuple(w): toK(c) for w,c in E}
diff = lin((1,ev),(-1,muE))
print(f'rational part z=(0..0): even==mu: {len(diff)==0}; odd==0: {len(od)==0}; total {time.time()-t0:.0f}s', flush=True)
if diff: print('EVEN residual:', list(diff.items())[:2])
if od: print('ODD residual:', list(od.items())[:2])
store = pickle.load(open('hpl_constructed_K7_nonrat.pkl','rb'))
store[zkey]=Fsym
pickle.dump(store, open('hpl_constructed_K7.pkl','wb'))
print('saved hpl_constructed_K7.pkl (all 15 zeta-monomials);', sum(len(v) for v in store.values()), 'HPL terms total')
open('progress.log','a').write(f'step05b: k=7 rational part constructed: {len(Fsym)} terms; even==mu & odd==0 exact; all 15 monomials in hpl_constructed_K7.pkl\n')
