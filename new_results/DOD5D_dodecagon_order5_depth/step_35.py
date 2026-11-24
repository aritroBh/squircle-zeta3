import pickle, sympy as sp, sys
sys.path.insert(0, "<project>/research/08_theory/E5c_nonlocal_z5")
import mpmath as mp
mp.mp.dps = 50
import symeng
# rebuild E5full
_t_, _al, _be, _Ls, _hs, _s0s = sp.symbols('t alpha beta L h s0', real=True)
_l0s, _l1s = sp.symbols('l0 l1')
_taus = sp.symbols('tau1:6'); _gams = sp.symbols('gam1:6')
_loc = {'t':_t_,'alpha':_al,'beta':_be,'L':_Ls,'h':_hs,'s0':_s0s,'l0':_l0s,'l1':_l1s}
for _i,_x in enumerate(_taus): _loc[f'tau{_i+1}']=_x
for _i,_x in enumerate(_gams): _loc[f'gam{_i+1}']=_x
_D = pickle.load(open("<project>/research/08_theory/E5_edge_law/steps/step_02_full.pkl",'rb'))
E5full = sp.sympify(_D['E5'], locals=_loc)

r3 = mp.sqrt(3)
M = []
for j in range(6):
    c, s = mp.cos(j*mp.pi/6), mp.sin(j*mp.pi/6)
    sg = 1 if j <= 3 else -1
    M.append((sg*c, sg*s))
y0 = (2-r3, mp.mpf(0)); y1 = ((2*r3-3)/2, (2-r3)/2)
P = symeng.edge_params(M, y0, y1)
t = mp.mpf('0.5')
# parameter values exactly as E5_sym computes them
c0, c1 = P['c0'], P['c1']
a = []; ap = []
for j in range(6):
    lin = c0[j]*(1-t)+c1[j]*t
    sg = 1 if (c0[j]+c1[j])>0 else -1
    a.append(abs(lin)); ap.append(sg*(c1[j]-c0[j]))
logs = [mp.log(a[j]) for j in range(6)]
l0v, l1v = logs[P['j0']], logs[P['j1']]
tauk = [sum(a[j]*logs[j]**k for j in P['nonloc']) for k in range(1,6)]
gamk = [sum(ap[j]*logs[j]**k for j in P['nonloc']) for k in range(1,6)]
vals = {**{_t_: t, _al: P['alpha'], _be: P['beta'], _Ls: P['L'], _hs: P['h'], _s0s: P['s0'],
           _l0s: l0v, _l1s: l1v},
        **{_taus[i]: tauk[i] for i in range(5)}, **{_gams[i]: gamk[i] for i in range(5)}}
v_sympy = E5full.subs(vals)
v_mpf = symeng.E5_sym(t, M, P, conv='B')
print("E5full.subs:", mp.nstr(mp.mpf(str(v_sympy)), 30))
print("E5_sym      :", mp.nstr(v_mpf, 30))
# Now compare against reconstruction from my parsed terms
terms = pickle.load(open("E5_terms.pkl","rb"))
allv = [t, P['alpha'], P['beta'], P['L'], P['h'], P['s0'], l0v, l1v] + tauk + gamk
tot = mp.mpf(0)
for exp, rat in terms:
    m = mp.mpf(1)
    for i, e in enumerate(exp):
        if e: m *= allv[i]**e
    tot += mp.mpf(str(rat))*m
print("my terms    :", mp.nstr(tot, 30))
