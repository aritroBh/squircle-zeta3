import sympy as sp, pickle, os
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd)
r5 = sp.sqrt(5)
rho = sp.sqrt(10-2*r5)
s1atom = sp.sqrt(sp.Rational(5,8) - r5/8)
s2atom = sp.sqrt(sp.Rational(5,8) + r5/8)
def canon(e):
    return sp.expand(e.xreplace({s1atom: rho/4, s2atom: rho*(1+r5)/8}))
def reduce_nf(e):
    """e in Q(sqrt5,rho) -> (a,b,c,d) rationals, e = a+b*r5+c*rho+d*rho*r5; works on expanded forms."""
    e = canon(e)
    # put rho as a true atom by replacing (10-2r5)^(1/2) occurrences: use Pow->symbol
    RH = sp.Symbol('RH'); R5 = sp.Symbol('R5')
    e2 = e.xreplace({rho: RH, r5: R5})
    a=b=c=d=sp.Rational(0)
    for tm in sp.Add.make_args(sp.expand(e2)):
        dd = tm.as_powers_dict()
        ir = int(dd.get(RH,0)); i5 = int(dd.get(R5,0))
        rest = tm/(RH**ir * R5**i5) if (ir or i5) else tm
        rat = sp.Rational(sp.simplify(rest))
        # reduce val = rat * RHO^ir * R5^i5 back into basis
        val = rat
        for _ in range(ir):
            val *= RH
        for _ in range(i5):
            val *= R5
        # now reduce: R5^2->5, RH^2->10-2R5
        val = sp.expand(val)
        changed = True
        while changed:
            changed = False
            val = sp.expand(val)
            val2 = val.xreplace({RH**2: 10-2*R5})
            if val2 != val: changed=True; val=val2
            val2 = val.xreplace({R5**2: 5})
            if val2 != val: changed=True; val=val2
        for tm2 in sp.Add.make_args(sp.expand(val)):
            dd2 = tm2.as_powers_dict()
            jr = int(dd2.get(RH,0)); j5 = int(dd2.get(R5,0))
            base = sp.Rational(tm2/(RH**jr*R5**j5) if (jr or j5) else tm2)
            if   (jr,j5)==(0,0): a+=base
            elif (jr,j5)==(0,1): b+=base
            elif (jr,j5)==(1,0): c+=base
            elif (jr,j5)==(1,1): d+=base
            else: raise ValueError(f"bad nf term {tm2}")
    return a,b,c,d
def is_zero(e): return reduce_nf(e)==(0,0,0,0)

phi = (1+r5)/2
n = [ ( (1+r5)/4, rho/4 ), ( (r5-1)/4, rho*(1+r5)/8 ), ( (1-r5)/4, rho*(1+r5)/8 ),
      ( -(1+r5)/4, rho/4 ), ( sp.Integer(-1), sp.Integer(0) ) ]
h  = (r5-1)/4            # apothem candidate = 1/(2*phi)
y0 = (r5-2)/rho          # half edge length candidate
V0 = (-h, y0); V1 = (-h, -y0)
t = sp.Symbol('t')
y  = (V0[0] + t*(V1[0]-V0[0]), V0[1] + t*(V1[1]-V0[1]))
print("kink checks (should be 0):")
print(" <n2,V0> =", reduce_nf(V0[0]*n[1][0]+V0[1]*n[1][1]))
print(" <n3,V1> =", reduce_nf(V1[0]*n[2][0]+V1[1]*n[2][1]))
# signs at t=1/2 (numeric) and affine forms
import mpmath as mp
sig = []
avals = []
for i in range(5):
    ex = y[0]*n[i][0] + y[1]*n[i][1]
    sg = 1 if float(sp.N(ex.subs(t, sp.Rational(1,2)))) > 0 else -1
    sig.append(sg)
    a = sp.expand(sg*ex)          # a_i(t) = |<n_i, y(t)>|
    avals.append(a)
    nf0 = reduce_nf(a.subs(t,0)); nf1 = reduce_nf(a.subs(t,1))
    print(f"a_{i+1}(t): nf(a(0))={nf0} nf(a(1))={nf1}  a(t)={sp.expand(a)}")
# sum = 1 (edge in boundary)
S = sp.expand(sum(avals))
print("sum a_i - 1 nf:", reduce_nf(S-1))
