
part2 = r'''

# ======================================================================
# Part 2: the octagon (exact)
# ======================================================================
# Families on 0<t<1:
#   G(t) = Li_21(-t, -(1-t)/(2t))        curve (x,xy) = (-t, (1-t)/2)
#   H(t) = Li_21((1-t)/2, -2t/(1-t))     curve (x,xy) = ((1-t)/2, -t)  [x<->p swap of G]
# Chain rule (Part 0) + Li_11 closed form give
#   G' = JG/t + Dm/(t(1+t)) + 2 Db/(1-t^2),   JG = (q-R) B + Dc - D1
#   H' = -JH/(1-t) - Dm/(t(1+t)) - 2 Db/(1-t^2),  JH = -lc q + Du - Dt
JG = (q-R)*B + Dc - D1
JH = -lc*q + Du - Dt
dG = JG/t + Dm/(t*(1+t)) + 2*Db/(1-t**2)
dH = -JH/(1-t) - Dm/(t*(1+t)) - 2*Db/(1-t**2)
# Closed forms (constants included; derivatives of constants vanish):
G_cl = 3*Z3/4 + Tt + Tm + Tmc - Tu - q*(Dt + Db/2 - Dh/2)
H_cl = T1 + Tb + Tw - Tc - sp.Symbol('T12') - B*(Dc + Db) + B**2*U/2
# T12 = Li3(1/2), a constant.
DERIV[sp.Symbol('T12')] = 0
# Function identities of Part 1 as substitution rules (proved there):
BASIC = {D1: Z - Dt - R*A,                 # E
         Dh: Z - Db - (A-C)*q,             # E at h=(1+t)/2  (1-h=b)
         Dmc: -Db - q**2/2,                # L2
         Du: Z - Dc - lc*U,                # E at u=2t/(1+t) (1-u=c)
         Dw: -Dm - B**2/2}                 # L1
OC = Dt - Dm + Dc + Db + q**2/2 - 3*Z/2 + R*lc   # = OC with Dmc eliminated via L2

def part2():
    print("=" * 74)
    print("PART 2: octagon derivative identities (exact)")
    rH = sp.factor((diff(H_cl) - dH).subs(BASIC))
    print("  d/dt H_closed - H'  (mod E,L1,L2) =", rH)
    assert rH == 0
    rG = sp.factor((diff(G_cl) - dG).subs(BASIC))
    print("  d/dt G_closed - G'  (mod E,L1,L2) =", sp.factor(rG), " = -OC/(1+t):",
          sp.factor(rG - (-OC/(1+t))))
    assert sp.factor(rG - (-OC/(1+t))) == 0
    print("  and OC = 0 by Part 1 (derivative 0 + pin at t->1-).")
    print("  constants:")
    print("   G(t)->0 as t->1- (|G|<=|y| zeta(2)/(1-|y|), y=-(1-t)/(2t)->0);")
    print("   G_closed(1) = 3Z3/4 + Li3(1) + Li3(-1) + Li3(0) - Li3(1) - 0")
    print("              = 3Z3/4 + Z3 - 3Z3/4 - Z3 = 0  [Li3(-1)=-3Z3/4, Part 1 DUP3]")
    print("   H(t)->0 as t->0+ (same estimate, y=-2t/(1-t)->0);")
    print("   H_closed(0) = Li3(1)+Li3(1/2)+Li3(0)-Li3(1)-Li3(1/2)-0+0 = 0 (exact)")
    print("  Hence G=G_closed and H=H_closed on (0,1).")
    print("PART 2 OK")


# ---------------- octagon specialisation at t = k = sqrt(2)-1 ----------------
def part2_specialisation():
    print("=" * 74)
    print("PART 2b: octagon specialisation at t = sqrt(2)-1 (exact over Q(sqrt(2)))")
    s = sp.sqrt(2)
    k, h, x2, d, x1 = s-1, s/2, 1-s/2, 2-s, 1-s
    args = {
        "b(k)=x2": (1-k)/2 - x2, "h(k)=sqrt(2)/2": (1+k)/2 - h,
        "c(k)=k": (1-k)/(1+k) - k, "u(k)=d": 2*k/(1+k) - d,
        "w(k)=x2": k/(1+k) - x2,
        "G second arg = -sqrt(2)/2": -(1-k)/(2*k) + h,
        "H second arg = -sqrt(2)": -2*k/(1-k) + s,
        "h(k)^2 = 1/2": h**2 - sp.Rational(1, 2),
        "k(1+sqrt2) = 1": k*(1+s) - 1,
        "1-k = sqrt(2) k": 1-k - s*k,
        "x2 = x1/(x1-1)": x2 - x1/(x1-1),
        "-sqrt2/2 = x2/x1": -s/2 - x2/x1,
        "-sqrt2 = x1/x2": -s - x1/x2,
    }
    for name, r in args.items():
        r = sp.simplify(r)
        print("   ", name, ": residual =", r)
        assert r == 0
    # silver dilog rules as specialisations of E and Landen (proved function ids):
    l, L = sp.log(2), sp.log(1+s)
    Z2s = sp.pi**2/6
    print("   silver rules (exact specialisations of E and Landen at points of (0,1)):")
    print("    D(d):  E at d=2-sqrt2: 1-d=k, log d = l/2 - L, log k = -L")
    print("    D(x2): Landen at x1: -x1/(1-x1) = x2, log(1-x1) = l/2")
    print("    D(h):  E at h=sqrt2/2: 1-h=x2, log h = -l/2, log x2 = -l/2 - L")
    print("    D(-h): Landen at k: -k/(1-k) = -h, log(1-k) = l/2 - L")
    checks = {"-x1/(1-x1) = x2": sp.simplify(-x1/(1-x1) - x2),
              "-k/(1-k) = -h": sp.simplify(-k/(1-k) + h),
              "log d = l/2 - L": sp.simplify(sp.log(d) - (l/2 - L)),
              "log k = -L": sp.simplify(sp.log(k) + L),
              "log x2 = -l/2 - L": sp.simplify(sp.log(x2) - (-l/2 - L)),
              "log(1-x1) = l/2": sp.simplify(sp.log(1-x1) - l/2)}
    for name, r in checks.items():
        print("     ", name, ": residual =", r)
        assert r == 0
    # compact specialised forms and exact comparison with OCT6 reductions.pkl
    reductions = pickle.loads((SRC_OCT / "reductions.pkl").read_bytes())
    l2, lu, z2, z3 = sp.symbols("l2 lu z2 z3")
    Dd, Dk, Dx1, Dx2, Dh, Dmh = sp.symbols("L2d L2k L2x1 L2x2 L2r2h L2mr2h")
    Tk, Tx1, Tx2, Td, Th, Tmh = sp.symbols("L3k L3x1 L3x2 L3d L3r2h L3mr2h")
    compact3 = 3*z3/4 + Tk + 2*Tx1 - Td + l2*Dk/2 + l2*(Dx2-Dh)/4
    compact4 = Td + 2*Tx2 - Tk - 4*Th - 4*Tmh - l2*(Dk+Dx2)/2 + l2**2*(l2/2-lu)/8
    rules = {Dd: z2-Dk+l2*lu/2-lu**2, Dx2: -Dx1-l2**2/8,
             Dh: z2+Dx1-l2*lu/2-l2**2/8, Dmh: -Dk-lu**2/2+l2*lu/2-l2**2/8}
    for key, compact in (("Lam", compact3), ("swap", compact4)):
        r = sp.expand((compact - reductions[key]).xreplace(rules))
        print("   ", key, ": compact specialisation - stored reduction (exact) =", r)
        assert r == 0
    print("   [compact forms: G_closed(k), H_closed(k) rewritten by the argument")
    print("    identities above; -log h(k) = l2/2, log(1+k) = l2/2,")
    print("    log u(k) = log d = l2/2 - L, and Li3(1/2) = 4Li3(h)+4Li3(-h) (DUP3).]")
    print("PART 2b OK")
'''
with open('verify_full.py','a') as f: f.write(part2)
print("part2 appended")
