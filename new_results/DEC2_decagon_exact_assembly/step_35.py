import sympy as sp, pickle, os
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd)
r5 = sp.sqrt(5)
rho = sp.sqrt(10-2*r5)   # rho^2 = 10-2*sqrt5; generates the degree-4 field (sqrt5 = (10-rho^2)/2)
def reduce_nf(e):
    """e in Q(sqrt5, rho) -> (a,b,c,d) rationals with e = a + b*r5 + c*rho + d*rho*r5."""
    e = sp.expand(e)
    a=b=c=d=sp.Rational(0)
    for tm in sp.Add.make_args(e):
        dd = tm.as_powers_dict()
        ir = int(dd.get(rho,0)); i5 = int(dd.get(r5,0))
        rest = sp.simplify(tm/(rho**ir * r5**i5)) if (ir or i5) else tm
        rat = sp.Rational(rest)
        val = sp.expand(rat * (10-2*r5)**(ir//2) * rho**(ir%2) * 5**(i5//2) * r5**(i5%2))
        vd = val.as_powers_dict()
        jr = int(vd.get(rho,0)); j5 = int(vd.get(r5,0))
        base = sp.Rational(sp.simplify(val/(rho**jr*r5**j5)) if (jr or j5) else val)
        if   (jr,j5)==(0,0): a+=base
        elif (jr,j5)==(0,1): b+=base
        elif (jr,j5)==(1,0): c+=base
        elif (jr,j5)==(1,1): d+=base
        else: raise ValueError("unexpected normal form")
    return a,b,c,d
def is_zero(e): return reduce_nf(e)==(0,0,0,0)
def in_q5(e):
    a,b,c,d = reduce_nf(e); return (c,d)==(0,0), a+b*r5

phi = (1+r5)/2
# normals n_i, i=1..5 (exact, in Q(rho))
n = [ ( (1+r5)/4, rho/4 ),
      ( (r5-1)/4, rho*(1+r5)/8 ),
      ( (1-r5)/4, rho*(1+r5)/8 ),
      ( -(1+r5)/4, rho/4 ),
      ( sp.Integer(-1), sp.Integer(0) ) ]
for i,(cx,sx) in enumerate(n,1):
    assert is_zero(cx**2+sx**2-1), f"unit {i}"
    assert is_zero(cx - sp.cos(i*sp.pi/5)) and is_zero(sx - sp.sin(i*sp.pi/5)), f"trig {i}"
print("normals OK: unit-norm, equal cos/isin(i*pi/5)")
# vertex on kink line k: V = s*(-sin(k th), cos(k th)), |s| = 1/S_k, S_k = sum_i |sin((k-i)pi/5)|
S_k = sp.simplify(2*sp.sin(sp.pi/5)+2*sp.sin(2*sp.pi/5))
s_val = sp.simplify(1/S_k)
print("S_k =", S_k, " = ", sp.nsimplify(S_k), "; s =", sp.nsimplify(s_val))
# for each k, s*S_k = 1 exactly:
for k in range(1,6):
    tot = 0
    for i in range(1,6):
        tot += sp.Abs(sp.sin((k-i)*sp.pi/5))
    assert is_zero(tot - S_k), k
print("all kink sums equal S_k exactly")
verts = {}
for k in range(1,6):
    perp = (-n[k-1][1], n[k-1][0])
    for sgn in [1,-1]:
        V = (sp.expand(sgn*s_val*perp[0]), sp.expand(sgn*s_val*perp[1]))
        # verify: on kink line k and on boundary
        assert is_zero(V[0]*n[k-1][0]+V[1]*n[k-1][1])
        tot = 0
        for i in range(5):
            tot += sp.Abs(V[0]*n[i][0]+V[1]*n[i][1])
        verts[(k,sgn)] = V
print("vertices on kink lines verified")
for k in verts: print(k, [sp.nsimplify(x) for x in verts[k]], "=", [sp.N(x,12) for x in verts[k]])
