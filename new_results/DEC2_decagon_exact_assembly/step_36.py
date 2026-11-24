import sympy as sp, pickle, os
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd)
r5 = sp.sqrt(5)
rho = sp.sqrt(10-2*r5)
s1atom = sp.sqrt(sp.Rational(5,8) - r5/8)   # sympy's sin(pi/5)
s2atom = sp.sqrt(sp.Rational(5,8) + r5/8)   # sympy's sin(2pi/5)
def canon(e):
    e = e.xreplace({s1atom: rho/4, s2atom: rho*(1+r5)/8})
    return sp.expand(e)
def reduce_nf(e):
    e = canon(e)
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
phi = (1+r5)/2
n = [ ( (1+r5)/4, rho/4 ),
      ( (r5-1)/4, rho*(1+r5)/8 ),
      ( (1-r5)/4, rho*(1+r5)/8 ),
      ( -(1+r5)/4, rho/4 ),
      ( sp.Integer(-1), sp.Integer(0) ) ]
for i,(cx,sx) in enumerate(n,1):
    assert is_zero(cx**2+sx**2-1), f"unit {i}"
    assert is_zero(cx - sp.cos(i*sp.pi/5)) and is_zero(sx - sp.sin(i*sp.pi/5)), f"trig {i}"
print("normals OK: unit norm, exact trig match")
# kink sums
S_k = 2*(rho/4) + 2*(rho*(1+r5)/8)   # 2 sin(pi/5)+2 sin(2pi/5)
for k in range(1,6):
    tot = 0
    for i in range(1,6):
        tot += sp.Abs(sp.sin((k-i)*sp.pi/5))
    assert is_zero(tot - S_k), k
s_val = sp.expand(1/S_k)
print("S_k exact (as nf):", reduce_nf(S_k), " s_val =", sp.nsimplify(s_val))
verts = {}
for k in range(1,6):
    perp = (-n[k-1][1], n[k-1][0])
    for sgn in [1,-1]:
        V = (sp.expand(sgn*s_val*perp[0]), sp.expand(sgn*s_val*perp[1]))
        assert is_zero(V[0]*n[k-1][0]+V[1]*n[k-1][1]), "on kink"
        # boundary: sum_i |<n_i,V>| = 1  (compute signed sum with signs from numeric eval)
        import mpmath as mp
        signs = [1 if float(sp.N(V[0]*n[i][0]+V[1]*n[i][1]))>0 else -1 for i in range(5)]
        tot = sum(signs[i]*(V[0]*n[i][0]+V[1]*n[i][1]) for i in range(5))
        assert is_zero(tot - 1), f"boundary {k},{sgn}"
        verts[(k,sgn)] = V
print("vertices verified: on kink line, on boundary sum=1")
R = sp.sqrt(verts[(2,-1)][0]**2 + verts[(2,-1)][1]**2)
for k in verts:
    print(k, [sp.nsimplify(x) for x in verts[k]])
