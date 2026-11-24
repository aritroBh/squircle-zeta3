# Verdict: D reduces modulo both primes using proved frame identities; exact certificate in progress

Weight-two calibration succeeded. Weight-three frame rows increase the modular
rank of the proved base, but none of the supplied R relations individually lies
in the enlarged span. The weight-four target is a member modulo both primes,
even after removing all R1–R5 lifts. No proof of Conjecture 4 is claimed until
the exact certificate passes an exact rational check.

| weight | base rank | + frame rows | quotient | target rank |
|---|---:|---:|---:|---:|
| 2 | 54 | 55 | 9 | 55 |
| 3 | 544 | 545 | 31 | 549 |
| 4, T3 (includes numerical R lifts) | 5074 | 5082 | 102 | 5103 |
| 4, T3noR (proved families) | 5059 | 5070 | 114 | 5103 |

Both rows agree at primes `1000000009` and `1000000033`, with square roots
`291087696` and `12641440`. R1–R5 are all nonmembers at both primes.

The direct series/G checks covered 25 identities at 65 decimal digits with maximum
residual `1.0385777e-66`. All 44 distinct weight-two rows and all 92 distinct
weight-three rows passed numerical verification (maximum absolute residuals
`1.2628346e-101` and `2.5924254e-101`). No row was discarded.

All 184 distinct weight-four frame rows passed with maximum absolute residual
`4.2057497e-101`. With all lower-weight frame lifts these give 3644 additional
rows. Every source row was checked against the supplied normal-form map at each
prime. Weight-four base ranks are reused from the supplied RREF artifacts; the
new ranks and target memberships are computed with FLINT.

For T3, the extra rows, extra rows plus P4, and extra rows plus P4 plus D have
normal-form ranks 8, 9, 9; for T3noR they are 11, 12, 12. The test uses the
combined D with its square-root coefficient embedded at each prime.

The simpler T0 system plus frame rows also contains D plus a multiple of P4
at the first prime (normal-form ranks 28, 29, 29). This permits a rational-only
certificate, avoiding both Kummer rows and R1–R5, if exact reconstruction succeeds.

The implementation uses exact rational arithmetic and the existing Cayley
conversion engine. Series and shuffle identities are checked independently with
mpmath direct sums. Progress and measured results are recorded in `progress.log`.
