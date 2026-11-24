import pickle, sympy as sp, time
_t_, _al, _be, _Ls, _hs, _s0s = sp.symbols('t alpha beta L h s0', real=True)
_l0s, _l1s = sp.symbols('l0 l1')
_taus = sp.symbols('tau1:6'); _gams = sp.symbols('gam1:6')
_loc = {'t':_t_,'alpha':_al,'beta':_be,'L':_Ls,'h':_hs,'s0':_s0s,'l0':_l0s,'l1':_l1s}
for _i,_x in enumerate(_taus): _loc[f'tau{_i+1}']=_x
for _i,_x in enumerate(_gams): _loc[f'gam{_i+1}']=_x
_D = pickle.load(open("<project>/research/08_theory/E5_edge_law/steps/step_02_full.pkl",'rb'))
E5full = sp.sympify(_D['E5'], locals=_loc)
allvars = [_t_,_al,_be,_Ls,_hs,_s0s,_l0s,_l1s]+list(_taus)+list(_gams)
t0 = time.time()
P5 = sp.Poly(E5full, *allvars)
terms = P5.terms()
print("terms:", len(terms), "parse time %.1fs" % (time.time()-t0))
# check weighted homogeneity: l0,l1 count 1; tau_k, gam_k count k
bad = []
degs = set()
for exp, cf in terms:
    w = exp[6] + exp[7] + sum((i+1)*exp[8+i] for i in range(5)) + sum((i+1)*exp[13+i] for i in range(5))
    degs.add(w)
    if w != 5: bad.append((exp, cf, w))
print("log-weight degrees present:", sorted(degs), " non-5 terms:", len(bad))
# t-degree and param-degree ranges
print("max t deg:", max(e[0] for e,_ in terms))
print("max alpha deg:", max(e[1] for e,_ in terms), "max beta:", max(e[2] for e,_ in terms))
print("L degs:", sorted(set(e[3] for e,_ in terms)), "h degs:", sorted(set(e[4] for e,_ in terms)), "s0 degs:", sorted(set(e[5] for e,_ in terms)))
print("sum of (L,h,s0) degs per term:", sorted(set(e[3]+e[4]+e[5] for e,_ in terms)))
print("alpha+beta degs:", sorted(set(e[1]+e[2] for e,_ in terms)))
pickle.dump(terms, open("E5_terms.pkl","wb"))
