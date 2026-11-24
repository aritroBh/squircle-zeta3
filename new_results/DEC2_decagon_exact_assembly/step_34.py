import sympy as sp, pickle, os
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd)
r5 = sp.sqrt(5)
rho = sp.sqrt(10-2*r5)   # degree-4 generator: rho^2 = 10-2*sqrt5 ; sqrt(10+2sqrt5)=4*sqrt5/rho
def reduce_nf(e):
    """normal form of e in Q(sqrt5, rho): returns (a,b,c,d) rationals with e = a + b*r5 + c*rho + d*rho*r5."""
    e = sp.expand(e)
    a=b=c=d=sp.Rational(0)
    for tm in sp.Add.make_args(e):
        dd = tm.as_powers_dict()
        ir = int(dd.get(rho,0)); i5 = int(dd.get(r5,0))
        rest = sp.simplify(tm/(rho**ir * r5**i5)) if (ir or i5) else tm
        rat = sp.Rational(rest)  # raises if not rational
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
# normals n_i = (cos(i pi/5), sin(i pi/5)), i=1..5
th = [i*sp.pi/5 for i in range(1,6)]
n = [(sp.cos(t), sp.sin(t)) for t in th]
print("normals exact:", n)
# vertices: V on kink line k means <n_k, V> = 0 and sum_i |<n_i,V>| = 1.
# V = s * (-sin(k th5), cos(k th5)); |s| = 1/sum_i |sin((k-i)th5)|
verts = {}
for k in range(1,6):
    nkperp = (-sp.sin(k*sp.pi/5), sp.cos(k*sp.pi/5))
    S = 0
    for i in range(1,6):
        S += sp.Abs(sp.sin((k-i)*sp.pi/5))
    # S is a sum of abs of sines -> all positive radicals; simplify exactly
    Sx = sp.simplify(2*sp.sin(sp.pi/5)+2*sp.sin(2*sp.pi/5))
    assert is_zero(S - Sx)
    s = sp.simplify(1/Sx)
    for sgn in [1,-1]:
        V = (sp.simplify(sgn*s*nkperp[0]), sp.simplify(sgn*s*nkperp[1]))
        verts[(k,sgn)] = V
# identify V0 (on kink line 2, upper-left) and V1 (on kink line 3, lower-left)
V0 = verts[(2,-1)]  # s negative? check components
print("verts[(2,-1)] =", V0)
print("verts[(2,1)]  =", verts[(2,1)])
print("verts[(3,-1)] =", verts[(3,-1)])
print("verts[(3,1)]  =", verts[(3,1)])
