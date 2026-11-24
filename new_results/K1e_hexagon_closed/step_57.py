
import mpmath as mp, json, pickle, sympy as sp, time
WORK='research/08_theory/K1e_hexagon_closed'
mp.mp.dps=110
t0=time.time()
direct = json.load(open(f"{WORK}/fresh_direct.json"))
X = mp.mpf(direct['1'])/mp.sqrt(13)
l3 = mp.log(3)
Li2m2 = mp.polylog(2,-2); Li3m2 = mp.polylog(3,-2); Li323 = mp.polylog(3, mp.mpf(2)/3)
z2 = mp.pi**2/6; z3 = mp.zeta(3)
basis = [Li3m2, Li323, Li2m2*l3, l3**3, Li2m2, l3**2, l3, z3, z2*l3, z2, mp.mpf(1)]
names = ['Li3(-2)','Li3(2/3)','Li2(-2)*log3','log^3 3','Li2(-2)','log^2 3','log3','zeta3','zeta2*log3','zeta2','1']
rel = mp.pslq([X]+basis, maxcoeff=10**12, maxsteps=10000)
print("PSLQ relation for E_arc/sqrt(13) (edge 1, r=3):")
if rel:
    print(rel, "  residual:", mp.nstr(abs(mp.fsum([mp.mpf(c)*b for c,b in zip(rel,[X]+basis)])),3))
    print("X =", " + ".join(f"({-mp.mpf(rel[i+1])/rel[0]})*{names[i]}" for i in range(11) if rel[i+1]!=0))
    coeff_pslq = {names[i]: sp.Rational(-rel[i+1], rel[0]) for i in range(11)}
    # compare with my formula coefficients for this edge
    clean = pickle.load(open(f"{WORK}/Earc_basis_clean.pkl","rb"))
    e = json.load(open(f"{WORK}/geometry_fresh.json"))[0]
    L2s=sp.Symbol('L2')
    subs={sp.Symbol('A1'):sp.Rational(e['A1']), sp.Symbol('B'):sp.Rational(e['B']), sp.Symbol('kap'):sp.Rational(e['kap']),
          L2s:sp.Integer(e['L2']), sp.Symbol('q0'):sp.Integer(e['q0']), sp.Symbol('V02'):sp.Integer(e['V02'])}
    # logs: lA0=log(A0/w0)=log((1/2)/3)=log(1/6); lA1=log((5/6)/5)=log(1/6); lB=log((1/6)/3)=log(1/18)
    la = {sp.Symbol('lA0'):-mp.log(6), sp.Symbol('lA1'):-mp.log(6), sp.Symbol('lB'):-mp.log(18)}
    mynames = {'Li3(x)':'Li3(-2)','Li3(x/(x-1))':'Li3(2/3)','Li2(x)*ln(1-x)':'Li2(-2)*log3','ln^3(1-x)':'log^3 3',
               'Li2(x)':'Li2(-2)','ln^2(1-x)':'log^2 3','ln(1-x)':'log3','zeta3':'zeta3','zeta2*ln(1-x)':'zeta2*log3','zeta2':'zeta2','1':'1'}
    # formula coefficient of atom = (num/L2^(5/2))/sqrt(13) with logs substituted; note log3 = -ln(1-x)? careful
    # basis atom ln(1-x) with x=-2 is ln(3); PSLQ basis uses log3 = ln 3. Same.
    maxdiff = sp.Integer(0)
    for name, ps in coeff_pslq.items():
        my = [n for n,v in mynames.items() if v==ps][0]
        num,den = clean[my]
        cv = sp.Rational(sp.N(num.subs(subs),30))/sp.Integer(13)**3   # /L2^(5/2)/sqrt(13) = /13^3
        # need exact: num has lA logs -> substitute numeric and PSLQ each? Instead evaluate coefficient of each log-monomial numerically
        print(f"{ps:>16}: pslq={ps}   (my coeff evaluated below)")
    # evaluate my coefficient numerically per atom with logs numeric and PSLQ against {log2,log3,1} combos is overkill;
    # direct numeric compare:
    maxd = mp.mpf(0)
    for name, psr in coeff_pslq.items():
        my = [n for n,v in mynames.items() if v==name][0]
        num,den = clean[my]
        cv = mp.mpf(sp.N(num.subs(subs), 60))
        cv = cv/13**3
        # substitute logs: num is polynomial in lA0,lA1,lB
        cv2 = mp.mpf(sp.N(sp.expand(num).subs({sp.Symbol('lA0'):-sp.log(6), sp.Symbol('lA1'):-sp.log(6), sp.Symbol('lB'):-sp.log(18)}).subs(subs), 60)) / 13**3
        d = abs(cv2 - mp.mpf(float(psr)))
        maxd = max(maxd, d)
        print(f"{name:>16}: pslq={str(psr):>12} mine={mp.nstr(cv2,12)} |d|={mp.nstr(d,3)}")
    print("max |pslq - formula coeff|:", mp.nstr(maxd,3))
print(f"({time.time()-t0:.0f}s)")
