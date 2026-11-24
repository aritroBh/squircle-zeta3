import pickle, mpmath as mp
mp.mp.dps=820
R=pickle.load(open('five_relations.pickle','rb')); F=pickle.load(open('free_vals_hi.pickle','rb'))

vals=[]
for w in R['freewords']:
    v=F['freevals'][tuple(w)]; vals.append(mp.mpc(mp.mpf(v[0]),mp.mpf(v[1])))
print('n words',len(vals),'sample',mp.nstr(vals[0],25))
for i,rel in enumerate(R['relations']):
    r=[int(c) for c in rel]; res=mp.fsum(c*v for c,v in zip(r,vals))
    print(f'R{i+1}: height {max(abs(c) for c in r)}, |residual| {mp.nstr(abs(res),3)}')
# rank of the 5 relations
import sympy as sp
M=sp.Matrix([[int(c) for c in rel] for rel in R['relations']]); print('rank of R1..R5 =', M.rank())
# independent spot check of one word at 30 digits: H(a,b,c) = int_{0<t3<t2<t1<1} f_a(t1) f_b(t2) f_c(t3), letter j: z^j/(1-z^j t), first index outermost
mp.mp.dps=30
z=mp.exp(2j*mp.pi/8); f=lambda j,t: z**j/(1-z**j*t)
w=R['freewords'][0]; print('word', w)
a,b,c=w
H1=lambda j,t: -mp.log(1-z**j*t)
inner=lambda t2: f(b,t2)*H1(c,t2)
H2=lambda t1: mp.quad(inner,[0,t1])
val=mp.quad(lambda t1: f(a,t1)*H2(t1),[0,1])
print('independent quad H', mp.nstr(val,20), ' stored', mp.nstr(vals[0],20), ' diff', mp.nstr(abs(val-vals[0]),3))
