"""CX1: exact extraction and independent checks; run with the project Python."""

from pathlib import Path
import hashlib
import pickle
import sys

import mpmath as mp
import sympy as sp


SOURCE = Path(__file__).resolve().parent.parent / "OCT6_octagon_closed"


def depth2(x, y, p=2, details=False):
    """Sum Li_{p,1}; tail <= rho**(N+1)/(1-rho), rho=max(|x|,|xy|)."""
    if p not in (1, 2):
        raise ValueError("Only Li_{1,1} and Li_{2,1} are used here")
    rho = max(abs(x), abs(x * y))
    if not 0 < rho < 1:
        if rho == 0:
            return mp.mpf(0)
        raise ValueError("The defining series requires |x|, |xy| < 1")
    target = mp.power(10, -mp.mp.dps - 5)
    nmax = int(mp.ceil(mp.log(target * (1 - rho)) / mp.log(rho)))
    # ponytail: an O(N) scaled recurrence avoids growing powers of y when |y| > 1.
    h, xy_power, total = mp.mpf(0), x * y, mp.mpf(0)
    for n in range(1, nmax + 1):
        total += h / n**p
        h = x * h + x * xy_power / n
        xy_power *= x * y
    return (total, nmax, rho**(nmax+1)/(1-rho)) if details else total


def G(t):
    return depth2(-t, -(1-t)/(2*t))


def H(t):
    return depth2((1-t)/2, -2*t/(1-t))


def G_closed(t):
    b, h = (1-t)/2, (1+t)/2
    c, u = b/h, 2*t/(1+t)
    return (3*mp.zeta(3)/4 + mp.polylog(3, t) + mp.polylog(3, -t)
            + mp.polylog(3, -c) - mp.polylog(3, u)
            - mp.log(h)*(mp.polylog(2, t)
                         + (mp.polylog(2, b)-mp.polylog(2, h))/2))


def H_closed(t):
    b, c, w = (1-t)/2, (1-t)/(1+t), t/(1+t)
    B, u = mp.log(1+t), 2*t/(1+t)
    return (mp.polylog(3, 1-t) + mp.polylog(3, b) + mp.polylog(3, w)
            - mp.polylog(3, c) - mp.polylog(3, mp.mpf(1)/2)
            - B*(mp.polylog(2, c)+mp.polylog(2, b))
            + B**2*mp.log(u)/2)


def named_values():
    s = mp.sqrt(2)
    args = {"x1": 1 - s, "x2": 1 - s / 2, "k": s - 1,
            "d": 2 - s, "r2h": s / 2, "mr2h": -s / 2}
    values = {"l2": mp.log(2), "lu": mp.log(1 + s),
              "z2": mp.zeta(2), "z3": mp.zeta(3)}
    values.update({f"L{n}{name}": mp.polylog(n, z)
                   for n in (2, 3) for name, z in args.items()})
    return values


def check_reductions():
    data = (SOURCE / "reductions.pkl").read_bytes()
    reductions = pickle.loads(data)
    mp.mp.dps = 220  # 200 requested digits plus guard digits.
    values = named_values()
    s = mp.sqrt(2)
    print("Exact extraction from", SOURCE / "reductions.pkl")
    print("SHA256 =", hashlib.sha256(data).hexdigest())
    for key, x, y in (("Lam", 1 - s, -s / 2), ("swap", 1 - s / 2, -s)):
        expression = reductions[key]
        print(key, "=", expression)
        rhs = sp.lambdify(sorted(expression.free_symbols, key=str), expression,
                          modules="mpmath")(*[values[str(v)] for v in
                          sorted(expression.free_symbols, key=str)])
        lhs, nmax, bound = depth2(x, y, details=True)
        residual = lhs - rhs
        print(key, "LHS =", mp.nstr(lhs, 205))
        print(key, "RHS =", mp.nstr(rhs, 205))
        print(key, "signed residual (220 dps) =", mp.nstr(residual, 16))
        print(key, "N =", nmax, "; absolute series tail bound =", mp.nstr(bound, 16))
        assert abs(residual) < mp.mpf("1e-200")
    for name, residual in (("G(k) - G_closed(k)", G(s-1)-G_closed(s-1)),
                           ("H(k) - H_closed(k)", H(s-1)-H_closed(s-1))):
        print(name, "(220 dps) =", mp.nstr(residual, 16))
        assert abs(residual) < mp.mpf("1e-200")
    D = lambda z: mp.polylog(2, z)
    l, L, Z = values["l2"], values["lu"], values["z2"]
    k, d, x1, x2, h = s-1, 2-s, 1-s, 1-s/2, s/2
    special = {
        "D(d) - (z2-D(k)+l*L/2-L^2)": D(d)-(Z-D(k)+l*L/2-L**2),
        "D(x2) - (-D(x1)-l^2/8)": D(x2)-(-D(x1)-l**2/8),
        "D(h) - (z2+D(x1)-l*L/2-l^2/8)": D(h)-(Z+D(x1)-l*L/2-l**2/8),
        "D(-h) - (-D(k)-L^2/2+l*L/2-l^2/8)": D(-h)-(-D(k)-L**2/2+l*L/2-l**2/8),
        "T(1/2) - 4*T(h) - 4*T(-h)": mp.polylog(3, mp.mpf(1)/2)-4*mp.polylog(3, h)-4*mp.polylog(3, -h),
    }
    for name, residual in special.items():
        print(name, "(220 dps) =", mp.nstr(residual, 16))
        assert abs(residual) < mp.mpf("1e-200")


def check_symbolic():
    """Differentiate the displayed formulas as rational polynomial algebra."""
    t = sp.symbols("t", positive=True)
    R, A, B, C, Z = sp.symbols("R A B C Z")
    Dt, Dm, Db, Dh, Dc, Dmc, D1, Du, Dw = sp.symbols("Dt Dm Db Dh Dc Dmc D1 Du Dw")
    Tt, Tm, Tmc, Tu, T1, Tb, Tw, Tc = sp.symbols("Tt Tm Tmc Tu T1 Tb Tw Tc")
    q, lc, U = B-C, A-B, C+R-B
    fg = Tt+Tm+Tmc-Tu-q*(Dt+Db/2-Dh/2)
    fh = T1+Tb+Tw-Tc-B*(Dc+Db)+B**2*U/2
    derivatives = {
        R:1/t, A:-1/(1-t), B:1/(1+t), Dt:-A/t, Dm:-B/t,
        D1:R/(1-t), Db:q/(1-t), Dh:-(A-C)/(1+t),
        Dc:2*U/(1-t**2), Dmc:-2*q/(1-t**2), Dw:B/(t*(1+t)),
        Tt:Dt/t, Tm:Dm/t, Tmc:-2*Dmc/(1-t**2), Tu:Du/(t*(1+t)),
        T1:-D1/(1-t), Tb:-Db/(1-t), Tw:Dw/(t*(1+t)), Tc:-2*Dc/(1-t**2),
    }
    diff = lambda e: sum(sp.diff(e, v)*dv for v, dv in derivatives.items())
    elementary = {
        "Euler derivative": diff(Dt+D1+R*A-Z),
        "Landen derivative at -t": diff(Dm+Dw+B**2/2),
        "Landen derivative at -c": diff(Dmc+Db+q**2/2),
        "odd complement derivative": diff(Dt-Dm+Dc-Dmc+R*lc-3*Z/2),
    }
    print("Exact symbolic checks (residuals must equal 0)")
    for name, residual in elementary.items():
        residual = sp.factor(residual)
        print(name, "=", residual)
        assert residual == 0
    x, y = sp.symbols("x y")
    z = -y*(1-x)/(1-y)
    residual = sp.factor((1-y)*(1-z)-(1-x*y))
    print("Li11 logarithm argument identity =", residual)
    assert residual == 0
    Jg, Jh = (q-R)*B+Dc-D1, -lc*q+Du-Dt
    dg = Jg/t+Dm/(t*(1+t))+2*Db/(1-t**2)
    dh = -Jh/(1-t)-Dm/(t*(1+t))-2*Db/(1-t**2)
    basic = {D1:Z-Dt-R*A, Dh:Z-Db-(A-C)*q, Dmc:-Db-q**2/2,
             Du:Z-Dc-lc*U, Dw:-Dm-B**2/2}
    O = Dt-Dm+Dc+Db+q**2/2-3*Z/2+R*lc
    for name, f, d, expected in (("G", fg, dg, -O/(1+t)), ("H", fh, dh, sp.Integer(0))):
        delta = sp.factor((diff(f)-d).subs(basic))
        residual = sp.factor(delta-expected)
        print(name, "derivative difference minus predicted factor =", residual)
        assert residual == 0
        print(name, "derivative difference =", "-O/(1+t)" if name == "G" else "0")
        residual = sp.factor(delta.subs(Dc, sp.solve(O, Dc)[0]))
        print(name, "derivative residual after odd complement =", residual)
        assert residual == 0

    s = sp.sqrt(2)
    k, h, x2, d = s-1, s/2, 1-s/2, 2-s
    args = {
        "b(k)=x2": (1-k)/2-x2,
        "h(k)=sqrt(2)/2": (1+k)/2-h,
        "c(k)=k": (1-k)/(1+k)-k,
        "u(k)=d": 2*k/(1+k)-d,
        "w(k)=x2": k/(1+k)-x2,
        "G second argument=-sqrt(2)/2": -(1-k)/(2*k)+h,
        "H second argument=-sqrt(2)": -2*k/(1-k)+s,
        "h(k)^2=1/2": h**2-sp.Rational(1,2),
        "k*(1+sqrt(2))=1": k*(1+s)-1,
        "1-k=sqrt(2)*k": 1-k-s*k,
        "(1+h)/h=1+sqrt(2)": (1+h)/h-(1+s),
    }
    for name, residual in args.items():
        residual = sp.simplify(residual)
        print(name, ": residual =", residual)
        assert residual == 0
    reductions = pickle.loads((SOURCE / "reductions.pkl").read_bytes())
    l, L, Z2, Z3 = sp.symbols("l2 lu z2 z3")
    Dd, Dk, Dx1, Dx2, Dh, Dmh = sp.symbols("L2d L2k L2x1 L2x2 L2r2h L2mr2h")
    Td, Tk, Tx1, Tx2, Th, Tmh = sp.symbols("L3d L3k L3x1 L3x2 L3r2h L3mr2h")
    compact3 = 3*Z3/4+Tk+2*Tx1-Td+l*Dk/2+l*(Dx2-Dh)/4
    compact4 = Td+2*Tx2-Tk-4*Th-4*Tmh-l*(Dk+Dx2)/2+l**2*(l/2-L)/8
    rules = {Dd:Z2-Dk+l*L/2-L**2, Dx2:-Dx1-l**2/8,
             Dh:Z2+Dx1-l*L/2-l**2/8, Dmh:-Dk-L**2/2+l*L/2-l**2/8}
    for key, compact in (("Lam", compact3), ("swap", compact4)):
        print(key, "compact specialization =", compact)
        residual = sp.expand((compact-reductions[key]).xreplace(rules))
        print(key, "compact specialization minus stored reduction =", residual)
        assert residual == 0


def check_functional():
    mp.mp.dps = 100
    D = lambda z: mp.polylog(2, z)
    T = lambda z: mp.polylog(3, z)
    print("Numerical functional checks: 100 decimal digits; signed residuals")
    print("Required absolute residual threshold = 1e-85")
    for j in (1, 2, 3, 4):
        t = mp.mpf(j)/5
        b, h, c, u, w = (1-t)/2, (1+t)/2, (1-t)/(1+t), 2*t/(1+t), t/(1+t)
        R, A, B, C = mp.log(t), mp.log(1-t), mp.log(1+t), mp.log(2)
        q, lc = mp.log(h), mp.log(c)
        Jg = (q-R)*B+D(c)-D(1-t)
        Jh = -lc*q+D(u)-D(t)
        dg = Jg/t+D(-t)/(t*(1+t))+2*D(b)/(1-t**2)
        dh = -Jh/(1-t)-D(-t)/(t*(1+t))-2*D(b)/(1-t**2)
        checks = {
            "Euler": D(t)+D(1-t)+R*A-mp.zeta(2),
            "Landen(-t)": D(-t)+D(w)+B**2/2,
            "Landen(-c)": D(-c)+D(b)+q**2/2,
            "odd complement": D(t)-D(-t)+D(c)-D(-c)+R*lc-3*mp.zeta(2)/2,
            "trilog duplication": T(t)+T(-t)-T(t**2)/4,
            "Li11 G formula": depth2(-t, -(1-t)/(2*t), 1)-Jg,
            "Li11 H formula": depth2(b, -2*t/(1-t), 1)-Jh,
            "G chain derivative": mp.diff(G, t)-dg,
            "H chain derivative": mp.diff(H, t)-dh,
            "G_closed derivative": mp.diff(G_closed, t)-dg,
            "H_closed derivative": mp.diff(H_closed, t)-dh,
            "G-G_closed": G(t)-G_closed(t),
            "H-H_closed": H(t)-H_closed(t),
        }
        for label, x, y in (("G", -t, -(1-t)/(2*t)), ("H", b, -2*t/(1-t))):
            checks[label+" partial_x"] = mp.diff(lambda a: depth2(a, y), x)-depth2(x, y, 1)/x
            checks[label+" partial_y"] = mp.diff(lambda a: depth2(x, a), y)-(D(x)-D(x*y)/y)/(1-y)
            checks[label+" Li11 derivative"] = mp.diff(lambda a: depth2(a, y, 1), x)+mp.log(1-x*y)/(1-x)
        print("t =", mp.nstr(t, 2))
        for name, residual in checks.items():
            print(" ", name, "=", mp.nstr(residual, 16))
            assert abs(residual) < mp.mpf("1e-85"), (t, name, residual)
    endpoints = {
        "T(-1)+3*zeta(3)/4": T(-1)+3*mp.zeta(3)/4,
        "G constant pin at t->1": 3*mp.zeta(3)/4+T(1)+T(-1)-T(1),
        "H constant pin at t->0": T(1)+T(mp.mpf(1)/2)-T(1)-T(mp.mpf(1)/2),
    }
    for name, residual in endpoints.items():
        print(name, "=", mp.nstr(residual, 16))
        assert abs(residual) < mp.mpf("1e-85")


if __name__ == "__main__":
    print("Python", sys.version.split()[0], "; mpmath", mp.__version__, "; sympy", sp.__version__)
    check_reductions()
    check_symbolic()
    check_functional()
    print("ALL CHECKS PASSED")
