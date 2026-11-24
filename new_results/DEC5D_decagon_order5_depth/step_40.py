
import pickle, time
import sympy as sp
base = "<project>/research/08_theory/"
E5str = pickle.load(open(base+"E5_edge_law/steps/step_02_full.pkl","rb"))['E5']

t, al, be, L, h, s0 = sp.symbols('t alpha beta L h s0')
l0, l1 = sp.symbols('l0 l1')
taus = sp.symbols('tau1:6'); gams = sp.symbols('gam1:6')
loc = {'t':t,'alpha':al,'beta':be,'L':L,'h':h,'s0':s0,'l0':l0,'l1':l1}
for i,x in enumerate(taus): loc[f'tau{i+1}']=x
for i,x in enumerate(gams): loc[f'gam{i+1}']=x
t0=time.time()
E5 = sp.sympify(E5str, locals=loc)
print("parse %.1fs"%(time.time()-t0))

# stage 1: s0 = L*sh, h = L*sqrt(hh)  ->  E5/L free of L
sh, hh = sp.symbols('sh hh')
E5x = E5.subs({s0: L*sh, h: L*sp.sqrt(hh)})
E5x = sp.expand(E5x)
E5oL = sp.expand(E5x/L)
# check no L or sqrt(hh) remains
assert not E5oL.has(L), "L remains!"
print("E5/L free of L; terms:", len(sp.Add.make_args(E5oL)))
import pickle as pk
pk.dump(str(E5oL), open("E5_dec_stage1.pkl","wb"))
print("saved stage1")
