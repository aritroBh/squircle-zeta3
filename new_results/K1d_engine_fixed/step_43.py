
import json, sympy as sp, mpmath as mp
mp.mp.dps = 130
geo = json.load(open('geometry.json'))
ver = json.load(open('t3_b1_110.json'))
eng = json.load(open('engine_arcs_fixed.json'))

z2,z3,l2 = sp.symbols('z2 z3 l2')
z2v = mp.pi**2/6; z3v = mp.zeta(3); l2v = mp.log(2)
pi2l2 = mp.pi**2*l2v; l2_3 = l2v**3; pi2 = mp.pi**2; l2_2 = l2v**2

# basis for PSLQ: [X, z3, pi^2 log2, log2^3, pi^2, log2^2, log2, 1]
print("FRESH PSLQ TEST, b=1 (engine prediction vs independent PSLQ of P3_arc/P0)")
print("="*100)
for i, arc in enumerate(ver['arcs']):
    e = geo['1'][i]; L = mp.sqrt(mp.mpf(e['L2']))
    P3 = mp.mpf(arc['P3_arc']); P0 = mp.mpf(arc['P0'])
    X = P3/P0
    basis_num = [X, z3v, pi2l2, l2_3, pi2, l2_2, l2v, mp.mpf(1)]
    # PSLQ at ~105 digits (data is 110 digits)
    rel = mp.pslq(basis_num, tol=mp.mpf('1e-95'), maxcoeff=10**12, maxsteps=10000)
    if rel is None:
        print(f"edge {i+1}: PSLQ FAILED"); continue
    rel = list(rel)
    # residual at full precision
    resid = mp.fsum(mp.mpf(c)*b for c,b in zip(rel,basis_num))
    # engine prediction: exact relation from closed form
    key = f"b1_e{i+1}"
    Ec = sp.sympify(eng[key]['Earc_closed'], locals={'z2':z2,'z3':z3,'l2':l2})
    XL = sp.expand(Ec/sp.sqrt(int(e['L2'])))          # predicted X = E/L
    # express in basis {z3, pi^2 l2 (=6 z2 l2), l2^3, pi^2 (=6 z2), l2^2, l2, 1}
    cof = { 'z3': sp.expand(XL).coeff(z3),
            'pi2l2': sp.expand(XL).coeff(z2*l2)*sp.Rational(1,1), }
    XLz = sp.expand(XL)
    c_z3   = XLz.coeff(z3)
    c_z2l2 = XLz.coeff(z2*l2)     # coefficient of z2*l2 = (pi^2/6) l2 -> basis pi2l2 coeff = c_z2l2/6
    c_l23  = XLz.coeff(l2**3)
    c_z2   = sp.expand(XLz.subs({z3:0, l2:0})) .coeff(z2)
    c_l22  = sp.expand(XLz.subs({z3:0, z2:0})).coeff(l2**2)
    c_l2   = sp.expand(XLz.subs({z3:0, z2:0})).coeff(l2) - 2*c_l22*0 - 3*c_l23*0  # careful
    # do it cleanly by successive extraction
    rem = XLz
    c_z3 = rem.coeff(z3); rem = sp.expand(rem - c_z3*z3)
    c_z2l2 = sp.expand(rem).coeff(z2*l2); rem = sp.expand(rem - c_z2l2*z2*l2)
    c_l23 = sp.expand(rem).coeff(l2**3); rem = sp.expand(rem - c_l23*l2**3)
    c_z2 = sp.expand(rem).coeff(z2); rem = sp.expand(rem - c_z2*z2)
    c_l22 = sp.expand(rem).coeff(l2**2); rem = sp.expand(rem - c_l22*l2**2)
    c_l2 = sp.expand(rem).coeff(l2); rem = sp.expand(rem - c_l2*l2)
    c_1 = sp.expand(rem)
    assert sp.expand(rem - c_1) == 0
    # basis order [X, z3, pi2*l2, l2^3, pi2, l2^2, l2, 1]; pi2l2 = 6*z2*l2, pi2 = 6*z2
    rats = [c_z3, c_z2l2/6, c_l23, c_z2/6, c_l22, c_l2, c_1]
    den = sp.ilcm(*[sp.denom(q) for q in rats])
    exact_rel = [int(den)] + [int(-den*q) for q in rats]   # den*X - den*sum = 0
    g = 0
    for v in exact_rel: g = sp.gcd(g, v)
    exact_rel = [v//g for v in exact_rel]
    # sign-normalize both
    def norm(v):
        v = list(v)
        for x in v:
            if x != 0:
                return [y*(-1 if x<0 else 1) for y in v]
    reln = norm(rel); exactn = norm(exact_rel)
    # proportional?
    prop = all(a*exactn[0]==b*reln[0] for a,b in zip(reln,exactn))
    # z3-coefficient of X from PSLQ: rel[0]*X + rel[1]*z3 + ... =0 -> X_z3 = -rel[1]/rel[0]
    pslq_cz3 = sp.Rational(-rel[1], rel[0])
    pred_cz3 = sp.sympify(eng[key]['cz3'])/sp.sqrt(int(e['L2']))
    Ze = sp.sympify(e['Z'])/sp.sqrt(int(e['L2']))
    Ne = sp.sympify(eng[key]['Ne'])/sp.sqrt(int(e['L2']))
    print(f"--- b=1 edge {i+1}: r={e['r']} ---")
    print("  PSLQ relation      :", reln, " residual:", mp.nstr(abs(resid),3))
    print("  engine-exact rel   :", exactn)
    print("  proportional       :", prop)
    print("  PSLQ  cz3 (of X=E/L):", pslq_cz3)
    print("  engine cz3 (of E/L) :", sp.simplify(pred_cz3), " = Z_e + N_e =", sp.simplify(Ze), "+", sp.simplify(Ne),
          "=", sp.simplify(Ze+Ne), " match:", sp.simplify(pslq_cz3-(Ze+Ne))==0)
    print("  cz3 of E_arc (x L)  :", eng[key]['cz3'], "   N_e =", eng[key]['Ne'])
    print()
