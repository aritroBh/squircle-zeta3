#!/usr/bin/env python3
"""Reproduce CX6 constants and certified table validation (no new packages).

Run with <project>/.venv/bin/python.
The supplied 254.23 main-term theorem is an input, not recertified here.
"""
from pathlib import Path
from time import monotonic
import mpmath as mp
from flint import acb, arb, ctx
from common import H

ROOT = Path(__file__).resolve().parents[3]
C_REST = 10_600_000
C_TOTAL = 1_210_000_000


def constants():
    p = arb.pi()
    a = p / (2 * arb(2).sqrt())
    x = 2 * (p / 4).exp()
    K = 64 / (p * a) * (a / 2 + p / 4).exp()
    T = 32 / (3 * p * a) * (3 * a / 2 + 3 * p / 4).exp()
    raw = 2 + T + K * ((p / 4).exp() + x.exp() + (arb(1).exp() * x).exp())
    delta = 2 - p.sqrt()
    conversion = (3 / (arb(1).exp() * delta)) ** 3
    total = arb("254.23") + C_REST * conversion
    for name, value in (
        ("pi", p), ("sqrt(pi)", p.sqrt()), ("a=pi/(2 sqrt(2))", a),
        ("x=2 exp(pi/4)", x), ("K", K), ("T", T),
        ("C_rest_raw", raw), ("C_rest - C_rest_raw", C_REST - raw),
        ("2-sqrt(pi)", delta), ("conversion_sup", conversion),
        ("C_total_raw_using_rounded_C_rest", total),
        ("C_total - C_total_raw", C_TOTAL - total),
    ):
        print(f"{name} = {value.str(30)}", flush=True)
        assert value.is_finite()
    assert delta > 0 and raw < C_REST and total < C_TOTAL
    print(f"CERTIFIED: C_rest={C_REST}, a_power=0, c=2; C_total={C_TOTAL}", flush=True)


def main_mp(k):
    """The requested integral, rotated by Cauchy's theorem; print 60 digits."""
    N = k + 1
    w = mp.expjpi(mp.mpf(1) / 4)
    rs = mp.sqrt(2 * N / mp.pi)

    def f(r):
        s = r * w
        return mp.expj(mp.pi * s / 2) * (s / (1 + s)) ** N * H(s) / s * w

    pts = [0] + [rs * t for t in (
        mp.mpf(1)/8, mp.mpf(1)/2, 1, 2, 4, 10)] + [mp.inf]
    return 2 / mp.pi * mp.quad(f, pts).imag


def main_arb(k):
    """Arb encloses the same integral; the omitted tail is proved in proof.md."""
    N = k + 1
    p = arb.pi()
    w = acb(1, 1) / arb(2).sqrt()

    def f(r, analytic):
        s = r * w
        # Meromorphic expression; the gamma pole at zero is cancelled exactly.
        return ((acb(0, 1) * p * s / 2).exp() * (s / (1 + s)) ** (N - 2) / (1 + s) ** 2
                * p.sqrt() / 2 * (s / 2 + 1).gamma()
                * (s / 2 + arb(3) / 2).rgamma() * w)

    # ponytail: fixed short cover for these three validations; use adaptive
    # endpoints only if the requested k-range changes.
    pts = [0, 1, 4, 8, 16, 32, 64, 128, 200]
    value = acb(0)
    for lo, hi in zip(pts, pts[1:]):
        piece = acb.integral(f, lo, hi, abs_tol=arb("1e-86"),
                             rel_tol=arb("1e-80"), eval_limit=100000)
        assert piece.is_finite(), (k, lo, hi, piece)
        value += piece
    a = p / (2 * arb(2).sqrt())
    R = arb(pts[-1])
    tail = 2 / p * arb(2).sqrt() * (-a * R).exp() / (a * R ** 2)
    # Use the certified upper endpoint as a symmetric outward error radius.
    enclosure = 2 / p * value.imag + arb(0, tail.upper())
    assert enclosure.is_finite() and enclosure.rad() < arb("1e-78")
    return enclosure, tail


def endpoint_audit(rows):
    """Print the elementary term which step_50 indexed incorrectly."""
    ell = arb(2).log()
    for k in (8, 10, 12):
        J = (-1) ** k * ell ** (k + 1) / arb(k + 1).fac()
        wrong = -arb(3) / 4 * (-ell) ** k / arb(k).fac()
        right = -arb(3) / 4 * J
        print(f"k={k} P_k(0)=J_k = {J}", flush=True)
        print(f"k={k} step_50 elementary term = {wrong}", flush=True)
        print(f"k={k} correct ray elementary term = {right}", flush=True)
        print(f"k={k} elementary discrepancy/|M_k| = {abs(wrong-right)/abs(rows[k])}", flush=True)


def representation_check(table_m0):
    """Check the new exact contour identity independently of the loose bound."""
    p = arb.pi()
    w = acb(1, 1) / arb(2).sqrt()

    def f(r, analytic):
        s = arb(1)/2 + r*w
        h = p.sqrt()/4 * (s/2).gamma() * (s/2+arb(3)/2).rgamma()
        # P_0(s)=integral_0^1 t^s/(1+t) dt, an exact digamma identity.
        pk = ((s/2+1).digamma() - ((s+1)/2).digamma())/2
        return h / (2*(p*s/2).cos()) * pk * w

    value = acb(0)
    pts = [0, 1, 4, 8, 16, 32, 64, 128, 200]
    for lo, hi in zip(pts, pts[1:]):
        piece = acb.integral(f, lo, hi, abs_tol=arb("1e-86"),
                             rel_tol=arb("1e-80"), eval_limit=100000)
        assert piece.is_finite(), (lo, hi, piece)
        value += piece
    a = p/(2*arb(2).sqrt())
    tail = 8/(p*a) * (-200*a).exp()  # |H|<=2, |K|<=2 exp(-a r), |P_0|<=1.
    result = -arb(2).log() + 2/p*value.imag + arb(0, tail.upper())
    difference = result-table_m0
    print(f"closed representation k=0 enclosure = {result}", flush=True)
    print(f"closed representation k=0 minus table = {difference}", flush=True)
    assert difference.contains(0) and difference.rad() < arb("1e-75")


def main():
    started = monotonic()
    ctx.dps = 90
    # Extra working digits preserve 60 significant digits after cancellation.
    mp.mp.dps = 90
    constants()
    rows = {}
    for line in (ROOT / "research/08_theory/D13f_sweep/mk_deg9.tsv").read_text().splitlines():
        fields = line.split("\t")
        if fields[0].isdigit():
            rows[int(fields[0])] = arb(fields[1])
    endpoint_audit(rows)
    representation_check(rows[0])
    for k in (100, 200, 400):
        approx = main_mp(k)
        certified, tail = main_arb(k)
        rest = abs(rows[k] - certified)
        N = arb(k + 1)
        bound = C_REST * (-2 * N.sqrt()).exp()
        print(f"k={k} supplied M_k = {rows[k]}", flush=True)
        print(f"k={k} M_main mpmath(90 working dps, 60 printed) = {mp.nstr(approx, 60)}", flush=True)
        print(f"k={k} M_main Arb enclosure = {certified}", flush=True)
        print(f"k={k} Arb main tail allowance = {tail}", flush=True)
        print(f"k={k} |M_k-M_main| Arb enclosure = {rest}", flush=True)
        print(f"k={k} rest bound = {bound}", flush=True)
        print(f"k={k} rest/bound = {rest/bound}", flush=True)
        print(f"k={k} validation summary: rest {rest.str(18)}; bound {bound.str(18)}", flush=True)
        assert rest < bound
        # Independent numerical agreement, not the source of the enclosure.
        assert abs(arb(mp.nstr(approx, 85)) - certified) < arb("1e-80")
    print(f"All checks passed; elapsed_seconds={monotonic()-started:.3f}", flush=True)


if __name__ == "__main__":
    main()
