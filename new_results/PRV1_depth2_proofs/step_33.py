
part3 = r'''

# ======================================================================
# Part 3: the decagon (exact)
# ======================================================================
# Families on 0<t<1 (same shape as the octagon's, with c=1):
#   G_dec(t) = Li_21(-t, -(1-t)/t)       curve (x,xy) = (-t, 1-t)
#   H_dec(t) = Li_21(1-t, -t/(1-t))      curve (x,xy) = (1-t, -t)  [x<->p swap]
Dv, Dw1, D1t2 = sp.symbols('Dv Dw1 D1t2')
T1t2, Tv, Tw1 = sp.symbols('T1t2 Tv Tw1')
DERIV.update({
    Dt2: -2*(A+B)/t,
    D1t2: 4*t*R/(1-t**2),
    Dv: B/(t*(1+t)),                    # v = t/(1+t)
    Dw1: (R-B)/(1+t),                   # w = 1/(1+t)
    Tt2: 2*Dt2/t,
    T1t2: -2*t*D1t2/(1-t**2),
    Tv: Dv/(t*(1+t)),
    Tw1: -Dw1/(1+t),
})
Li11_Gd = -R*B + D1 - D1t2
Li11_Hd = -R*A + Dt - Dt2
dGd = Li11_Gd/t + Dm/t + D1/(1-t)
dHd = -Li11_Hd/(1-t) - D1/(1-t) - Dm/t
G_dec_cl = Z3 - R*Dt2/2 + Tt2 - 2*Tt - T1
H_dec_cl = T1t2 - Tw1 + Tv - B*(Dm/2 + 3*D1/2 + D1t2/4 + Z/4)
BASIC_DEC = {Dt2: 2*Dt + 2*Dm,                       # DUP2
             D1: Z - Dt - R*A,                       # E
             D1t2: Z - (2*Dt+2*Dm) - 2*R*(A+B),      # E at t^2
             Dv: -Dm - B**2/2,                       # L1 (Landen at -t)
             Dw1: Z - (-Dm - B**2/2) + B*(R-B)}      # E at v (1-v = w)

def part3():
    print("=" * 74)
    print("PART 3: decagon derivative identities (exact)")
    rG = sp.factor((diff(G_dec_cl) - dGd).subs(BASIC_DEC))
    print("  d/dt G_dec_closed - G_dec' (mod DUP2,E,L1) =", rG)
    assert rG == 0
    rH = sp.factor((diff(H_dec_cl) - dHd).subs(BASIC_DEC))
    print("  d/dt H_dec_closed - H_dec' (mod DUP2,E,L1) =", rH)
    assert rH == 0
    print("  constants:")
    print("   G_dec(t)->0 as t->1- (|G|<=|y| zeta(2)/(1-|y|), y=-(1-t)/t->0);")
    print("   G_dec_closed(1) = Z3 - 0 + Li3(1) - 2 Li3(1) - Li3(0) = 0 (exact)")
    print("   H_dec(t)->0 as t->0+ (same estimate, y=-t/(1-t)->0);")
    print("   H_dec_closed(0) = Li3(1)-Li3(1)+Li3(0)-0 = 0 (exact)")
    print("  Hence G_dec=G_dec_closed and H_dec=H_dec_closed on (0,1).")
    print("PART 3 OK")


# ---- 3-term Landen trilog identity (function identity, exact) ----
Dml, Tml = sp.symbols('Dml Tml')     # Li2(-t/(1-t)), Li3(-t/(1-t))
DERIV[Dml] = A/(t*(1-t))
DERIV[Tml] = Dml/(t*(1-t))
T3 = Tt + T1 + Tml - (Z3 + Z*A - sp.Rational(1,2)*R*A**2 + A**3/6)

def part3_landen3():
    print("=" * 74)
    print("PART 3b: 3-term Landen trilog identity (exact)")
    r = sp.factor(diff(T3).subs({D1: Z - Dt - R*A, Dml: -Dt - A**2/2}))
    print("  d/dt[Li3(t)+Li3(1-t)+Li3(-t/(1-t)) - (Z3 + Z A - R A^2/2 + A^3/6)]")
    print("  after E and Landen =", r)
    assert r == 0
    print("  pin at t->0+: 0 + Z3 + 0 - (Z3 + 0 - 0 + 0) = 0 (exact); identity holds on (0,1).")
    print("PART 3b OK")


# ---- golden lattice and specialisation (exact over Q(sqrt(5))) ----
def part3_specialisation():
    print("=" * 74)
    print("PART 3c: golden lattice + decagon specialisation at t = 1/phi (exact)")
    s5 = sp.sqrt(5)
    phi = (1+s5)/2
    l, z2, z3 = sp.symbols('l z2 z3')
    L2a, L2b, L2c = sp.symbols('L2a L2b L2c')   # Li2(1/phi^2), Li2(1/phi), Li2(-1/phi)
    L3a, L3b, L3c = sp.symbols('L3a L3b L3c')   # Li3(1/phi^2), Li3(1/phi), Li3(-1/phi)
    # argument identities over Q(sqrt(5))
    argid = {
        "1 - 1/phi^2 = 1/phi": sp.simplify(1 - 1/phi**2 - 1/phi),
        "-(1/phi^2)/(1-1/phi^2) = -1/phi": sp.simplify(-(1/phi**2)/(1-1/phi**2) + 1/phi),
        "1/(1+1/phi) = 1/phi": sp.simplify(1/(1+1/phi) - 1/phi),
        "(1/phi)/(1+1/phi) = 1/phi^2": sp.simplify((1/phi)/(1+1/phi) - 1/phi**2),
        "log(1+1/phi) = log phi": sp.simplify(sp.log(1+1/phi) - sp.log(phi)),
        "(1/phi)^2 = 1 - 1/phi": sp.simplify((1/phi)**2 - (1-1/phi)),
    }
    for name, r in argid.items():
        print("   ", name, ": residual =", r)
        assert r == 0
    # golden dilogs: E @ 1/phi^2, Landen @ 1/phi^2, DUP2 @ 1/phi
    sold = sp.solve([L2a + L2b - (z2 - 2*l**2),
                     L2c + L2a + l**2/2,
                     L2a - 2*L2b - 2*L2c], [L2a, L2b, L2c], dict=True)[0]
    print("   golden dilogs:", sold)
    # golden trilogs: 3-term Landen @ 1/phi^2, DUP3 @ 1/phi
    solt = sp.solve([L3a + L3b + L3c - (z3 - l*z2 + 5*l**3/6),
                     L3a - 4*L3b - 4*L3c], [L3a, L3b, L3c], dict=True)[0]
    print("   Lewin: Li3(1/phi^2) =", solt[L3a])
    print("   sum relation: Li3(1/phi)+Li3(-1/phi) =", sp.expand(solt[L3b] + solt[L3c]))
    S = sp.expand(solt[L3b] + solt[L3c])
    # specialisation
    Lam = sp.expand((z3 + l*L2a/2 - 2*L3b).subs(sold))          # G_dec_closed(1/phi)
    print("   Lambda_dec = G_dec_closed(1/phi) =", Lam, "  [with Li3(1/phi) free]")
    Lam_DEC1 = 2*z3/5 + 4*l*z2/5 - l**3 + 3*L3c + L3b
    r = sp.simplify(sp.expand(Lam - Lam_DEC1).subs(L3b, S - L3c))
    print("   Lambda_dec - DEC1 form (mod sum relation) =", r)
    assert r == 0
    swap = sp.expand((L3a - l*(L2c/2 + 3*L2a/2 + L2b/4 + z2/4)).subs(sold))
    print("   swap_dec = H_dec_closed(1/phi) =", swap)
    swap_full = sp.expand(swap.subs(L3a, solt[L3a]))
    print("   swap_dec fully reduced (via Lewin) =", swap_full)
    swap_DEC1 = 3*z3/5 - 7*l*z2/5 + 2*l**3 + L3c + L3b
    r = sp.simplify(sp.expand(swap_full - swap_DEC1).subs(L3b, S - L3c))
    print("   swap_dec - DEC1 form (mod sum relation) =", r)
    assert r == 0
    print("PART 3c OK")
'''
with open('verify_full.py','a') as f: f.write(part3)
print("part3 appended")
