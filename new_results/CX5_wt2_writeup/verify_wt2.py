"""Exact checks accompanying proof.md; finite checks do not replace its proofs.

Run from any directory with the project .venv/bin/python -B and this file.
Reads the supplied trusted pickle and fastred.py; writes no files.
"""

import sys
sys.dont_write_bytecode = True

import pickle
from fractions import Fraction
from functools import lru_cache
from math import comb, factorial, prod
from pathlib import Path
from time import monotonic

import sympy as sp

SOURCE = Path(__file__).resolve().parent.parent / "WT2b_stabilization"
sys.path.insert(0, str(SOURCE))
from fastred import IQ, QSI, ZERO, pf_zero

TAU1 = QSI(Fraction(-1, 2), 0, Fraction(-1, 2))
TAU2 = TAU1.conj_i()
SQRT2 = QSI(0, 1)
failures = []
started = monotonic()


def check(label, run):
    try:
        detail = run()
        print(f"PASS: {label}" + (f"; {detail}" if detail else ""), flush=True)
    except Exception as exc:
        failures.append(label)
        print(f"FAIL: {label}: {type(exc).__name__}: {exc}", flush=True)


@lru_cache(None)
def K(d):
    odd = prod(range(1, 2 * d - 2, 2))
    return (IQ ** ((d + 1) % 4)) * QSI(Fraction((-1) ** d * odd, 2 * factorial(d))) / SQRT2


@lru_cache(None)
def rho(r, s):
    return K(r + s) * ((-1) ** r * (1 + (r > s)) * comb(r + s, r))


def hat(p):
    return rho(p + 1, p) * (p + 1) * TAU1 / (4 * p - 1)


with (SOURCE / "record_W2b.pkl").open("rb") as handle:
    record = pickle.load(handle)


def h(k, r, s):
    if r + s > k:
        return ZERO
    poles = record[k, (r, s)]["H"][1]
    order = 2 * (r + s) - 2
    return poles[order - 1] if len(poles) >= order else ZERO


def predicted(k, r, s, corrected=True):
    d = r + s
    feed = (r + 1) * TAU1 * h(k, r + 1, s)
    # The adjacent term already contains the ordinary second feed.
    if r >= s + 2 or (not corrected and r > s):
        feed += (s + 1) * TAU2 * h(k, r, s + 1)
    if s == r - 1:
        diagonal = h(k, r, r)
        feed += r * TAU2 * (diagonal - diagonal.conj_i())
    return (rho(r, s) if k == d else ZERO) + feed / (2 * d - 1)


def recorded_recursion():
    expected = {(k, (r, d - r)) for k in range(1, 13)
                for d in range(1, k + 1) for r in range((d + 1) // 2, d + 1)}
    assert {key for key in record if key[0] <= 12} == expected, "missing or unexpected keys"
    top_count = 0
    for k, (r, s) in sorted(expected):
        data = record[k, (r, s)]
        if r + s == 1:
            assert data["H"] == pf_zero() and data["Gm"] == pf_zero(), (k, r, s)
            continue  # Degree one is analytic; the top-pole recursion starts at d=2.
        assert h(k, r, s) == predicted(k, r, s), (k, r, s)
        assert h(k, r, s) == (rho(r, s) if r > s or k == r + s else hat(r)), (k, r, s)
        for poles in data["H"][1:3]:
            assert all(c.is_zero() for c in poles[2 * (r + s) - 2:]), (k, r, s)
        top_count += 1
    return f"all {top_count} keys with d>=2, including initial keys; 12 degree-one keys analytic"


def printed_failure():
    actual = h(4, 2, 1)
    old = predicted(4, 2, 1, corrected=False)
    assert actual == QSI(0, Fraction(-3, 4))
    assert actual - old == QSI(0, Fraction(3, 16), 0, Fraction(3, 16))
    assert actual != old
    return "k=4, key (2,1): actual - printed = 3*sqrt(2)*(1+i)/16 != 0"


def identity_1():
    for d in range(2, 201):
        assert (d + 1) * K(d + 1) == -IQ * (2 * d - 1) * K(d), d


def identity_2():
    count = 0
    for d in range(2, 201):
        for s in range((d - 2) // 2 + 1):
            r = d - s
            assert (2 * d - 1) * rho(r, s) == (
                (r + 1) * TAU1 * rho(r + 1, s) + (s + 1) * TAU2 * rho(r, s + 1)
            ), (d, r, s)
            count += 1
    return f"{count} nonadjacent keys"


def identity_3():
    for p in range(2, 201):
        diagonal = rho(p, p)
        assert (4 * p - 3) * rho(p, p - 1) == (
            (p + 1) * TAU1 * rho(p + 1, p - 1)
            + p * TAU2 * (diagonal - diagonal.conj_i())
        ), p


def identity_4():
    for p in range(1, 201):
        diagonal, stable = rho(p, p), hat(p)
        assert stable == (1 - IQ) * diagonal, p
        assert stable - stable.conj_i() == diagonal - diagonal.conj_i(), p


def finite_differences():
    for d in range(41):
        for m in range(d + 1):
            actual = sum((-1) ** j * comb(d, j) * j ** m * (d - j) ** (d - m)
                         for j in range(d + 1))
            assert actual == (-1) ** m * factorial(d), (d, m)
    return "all 861 pairs (d,m), including d=0"


v, eps = sp.symbols("v epsilon")
sqrt2 = sp.sqrt(2)
a, b = 1 + sqrt2, sqrt2 - 1
D1 = 1 / (v - b) + 1 / (v + a) - 1 / v
D2 = 1 / (v - a) + 1 / (v + b) - 1 / v
g20 = (1 + sqrt2) / 4 - sqrt2 / (2 * (v ** 2 + 1))
g11 = -sqrt2 * (v ** 2 + 2 * v - 1) / (4 * (v ** 2 + 1))
R10 = (1 - 1 / sqrt2) / v - 2 / (v + a)
Rw2 = 1 / (2 * v) + 1 / (2 * (v - a)) - 1 / (v + a) - 1 / (2 * (v + b))


def zero(expression):
    return sp.cancel(expression, extension=sqrt2) == 0


def epsilon_identity():
    G20 = eps ** 2 * g20 / (1 - eps)
    G11 = -eps ** 2 * g20 + eps ** 3 * g11 / (1 - eps)
    W10 = eps * R10 - 2 * G20 * D1 - G11 * D2 - (G11 * D1).subs(v, 1 / v) / v ** 2
    assert zero((1 - eps) * (W10 - eps * R10) - eps ** 2 * Rw2)
    # Also check that the primitives used here are the normalized double-pole primitives.
    for leading, primitive in [(rho(2, 0), g20), (rho(1, 1), -g20), (hat(1), g11)]:
        H = leading.sympy() / (v - sp.I) ** 2 + leading.conj_i().sympy() / (v + sp.I) ** 2
        assert zero(sp.diff(primitive, v) - H)
        assert sp.simplify(primitive.subs(v, b)) == 0
    return "symbolic in epsilon and v; normalized primitives checked too"


def integral_identity():
    # Exact analytic certificate, not numerical quadrature. See the last paragraph of proof.md.
    qv = (v ** 2 + 1) / (2 * sqrt2 * v)
    assert zero(Rw2 - D2 / 2 + sp.diff(qv, v) / (1 + qv))
    t, z = sp.symbols("t z", positive=True)
    qt = sp.sqrt(1 - 2 * t * (1 - t))
    f = sp.log((1 + qt) / 2)
    assert f.subs(t, 0) == f.subs(t, 1) == 0
    assert sp.limit(f / t, t, 0, dir="+") == -sp.Rational(1, 2)
    generating_function = sp.log((1 + sp.sqrt(1 - z)) / 2)
    derivative = (1 - 1 / sp.sqrt(1 - z)) / (2 * z)
    assert sp.simplify(sp.diff(generating_function, z) - derivative) == 0
    assert generating_function.subs(z, 0) == 0
    n = sp.symbols("n", integer=True, positive=True)
    # beta(n,n+1)=(n-1)! n!/(2n)!; the coefficient is 1/(2^(n+1)*n^2) for every n.
    coefficient = (sp.binomial(2 * n, n) * sp.factorial(n - 1) * sp.factorial(n)
                   / (2 * n * 2 ** n * sp.factorial(2 * n)))
    assert sp.simplify(coefficient - 1 / (2 ** (n + 1) * n ** 2)) == 0
    result = sp.pi ** 2 / 12 - sp.polylog(2, sp.Rational(1, 2)) / 2
    assert sp.simplify(result - sp.pi ** 2 / 24 - sp.log(2) ** 2 / 4) == 0
    return "pi^2/24 + log(2)^2/4, exact binomial/beta certificate for general n"


if __name__ == "__main__":
    check("corrected recursion and stabilization vs reducer data, k<=12", recorded_recursion)
    check("printed uncorrected recursion fails at the required key", printed_failure)
    check("I1 exactly, 2<=d<=200", identity_1)
    check("I2 exactly, 2<=d<=200", identity_2)
    check("I3 exactly, 2<=p<=200 (in particular all degrees <=200)", identity_3)
    check("I4 exactly, 1<=p<=200 (in particular all degrees <=200)", identity_4)
    check("finite-difference identity exactly, d<=40", finite_differences)
    check("(1-epsilon)(W10-epsilon*R10)=epsilon^2*Rw2", epsilon_identity)
    check("integral of Rw2*log(t)", integral_identity)
    elapsed = monotonic() - started
    print(f"{'FAIL' if failures else 'PASS'}: {9 - len(failures)}/9 checks, {elapsed:.2f}s", flush=True)
    sys.exit(bool(failures))
