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
t0=time.time()
E5exp = sp.expand(E5full)
print("expand: %.1fs" % (time.time()-t0))
terms = []
for term in sp.Add.make_args(E5exp):
    cf, mon = term.as_independent(*allvars, as_Mul=True)
    exp = [0]*18
    for f in sp.Mul.make_args(mon):
        if f.is_Symbol:
            exp[allvars.index(f)] += 1
        elif f.is_Pow and f.exp.is_integer:
            exp[allvars.index(f.base)] += int(f.exp)
        else:
            raise ValueError(f"unexpected factor {f}")
    terms.append((tuple(exp), sp.Rational(cf)))
print("terms:", len(terms), "parse time %.1fs" % (time.time()-t0))
degs = set(); ldegs = set(); pdegs = set()
for exp, cf in terms:
    w = exp[6] + exp[7] + sum((i+1)*exp[8+i] for i in range(5)) + sum((i+1)*exp[13+i] for i in range(5))
    degs.add(w); ldegs.add(exp[3]+exp[4]+exp[5]); pdegs.add(exp[1]+exp[2])
print("log-weight degs:", sorted(degs))
print("L+h+s0 degs:", sorted(ldegs), " alpha+beta degs:", sorted(pdegs))
print("L degs:", sorted(set(e[3] for e,_ in terms)), " max t deg:", max(e[0] for e,_ in terms))
pickle.dump(terms, open("E5_terms.pkl","wb"))
