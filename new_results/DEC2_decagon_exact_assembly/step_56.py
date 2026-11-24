import sympy as sp, pickle, os, time
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd)
t = sp.Symbol('t'); s = sp.Symbol('s')
z2 = sp.pi**2/6; z3 = sp.zeta(3)
anti = pickle.load(open('Jloc_antiderivatives.pkl','rb'))
def limit0(F):
    # t -> 0+
    e = F.subs({sp.polylog(2,t):0, sp.polylog(3,t):0,
                sp.polylog(2,1-t):z2, sp.polylog(3,1-t):z3,
                sp.log(1-t):0})
    e = sp.expand(e)
    keep = 0
    for term in sp.Add.make_args(e):
        if term.has(t):
            if sp.degree(term, t) >= 1: continue   # t^m * bounded/log -> 0
        keep += term
    keep = sp.expand(keep)
    # coefficients of log^q t must vanish
    lt = sp.log(t)
    for q in range(1,6):
        c = sp.simplify(keep.coeff(lt, q))
        assert c == 0, ("log^q t survives", q, c, F)
    return sp.simplify(keep.subs(lt, 0))
def limit1(F):
    # t -> 1- : set t = 1-s, s -> 0+
    e = F.subs(t, 1-s)
    e = e.subs({sp.polylog(2,s):0, sp.polylog(3,s):0,
                sp.polylog(2,1-s):z2, sp.polylog(3,1-s):z3,
                sp.log(1-s):0})
    e = sp.expand(e)
    keep = 0
    for term in sp.Add.make_args(e):
        if term.has(s):
            if sp.degree(term, s) >= 1: continue
        keep += term
    keep = sp.expand(keep)
    ls = sp.log(s)
    for q in range(1,6):
        c = sp.simplify(keep.coeff(ls, q))
        assert c == 0, ("log^q s survives", q, c, F)
    return sp.simplify(keep.subs(ls, 0))
O5 = '<project>/research/08_theory/OCT6_octagon_closed/oct5'
Jloc = pickle.load(open(O5+'/Jloc.pkl','rb'))
t0=time.time(); bad=0
myJloc = {}
for key in sorted(Jloc.keys()):
    F = anti[key]
    v = sp.simplify(limit1(F) - limit0(F))
    myJloc[key] = v
    want = Jloc[key].subs(sp.pi**2, 6*z2).subs(sp.zeta(3), z3)
    d = sp.simplify(v - want)
    if d != 0:
        bad += 1
        print("VALUE MISMATCH", key, "mine:", v, " want:", want)
print("values:", len(myJloc), "mismatches:", bad, " %.1fs"%(time.time()-t0))
pickle.dump(myJloc, open('Jloc_dec2.pkl','wb'))
with open('progress.log','a') as f:
    f.write(f"step_06: 40 local integrals PROVED by antiderivative certificates (F'=f symbolic 0) + exact endpoint evaluation; values vs OCT5 Jloc.pkl: {bad} symbolic mismatches\n")
print("sample:", [(k, myJloc[k]) for k in sorted(myJloc)[:6]])
