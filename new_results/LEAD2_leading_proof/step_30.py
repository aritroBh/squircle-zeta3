import os, time
os.chdir('<project>/research/08_theory/LEAD2_leading_proof')
p2 = open('proof_part2.md').read()

old55 = """    K_d/d! = (1/2)[z^d](pi/sin pi z) + (1/2)[z^d](pi/sin pi z)·(-1)·... = (1/4)[z^{d+1}] WW(z).

Precisely: `K_d/d! = T_d/d! - (u m_e part)/d! = (1/2)[z^d] pi/sin - (1/2)[z^d] pi/sin
+ (1/4)[z^{d+1}] WW(z) = (1/4)[z^{d+1}](E(z)+E(-z)) = (1/4)(e_{d+1} + (-1)^{d+1} e_{d+1}) = e_{d+1}/2`
since `d+1` is even. ∎"""
new55 = """    K_d/d! = T_d/d! - (u m_e part)/d!
           = (1/2)[z^d](pi/sin pi z) - (1/2)[z^d](pi/sin pi z) + (1/4)[z^{d+1}] WW(z)
           = (1/4)[z^{d+1}](E(z)+E(-z))
           = (1/4)(e_{d+1} + (-1)^{d+1} e_{d+1}) = e_{d+1}/2,

since `d+1` is even. ∎"""
assert old55 in p2
p2 = p2.replace(old55, new55)

old61 = p2[p2.find("**Lemma 6.1"):p2.find("**Theorem (the law)**")]
new61 = """**Lemma 6.1 (zigzag signs).** `delta_d := [x^d](sech x + tanh x) = eps_{d-1} A_d/d!` for all `d`,
and `eps_m A_{m+1}/m! = [z^m] (1-sinh z)/cosh^2 z`.

**Proof.** `(sec + tan)(ix) = sec(ix) + tan(ix) = sech x + i tanh x` (using `sec(ix) = sech x`,
`tan(ix) = i tanh x`). Comparing coefficients in `sum_d A_d (ix)^d/d! = sech x + i tanh x`:
`[x^d] sech x = Re(A_d i^d)/d! = A_d cos(pi d/2)/d!` (even `d`), `[x^d] tanh x = Im(A_d i^d)/d!
= A_d sin(pi d/2)/d!` (odd `d`). Since `sech` is even and `tanh` odd, `delta_d =
(-1)^{d/2} A_d/d!` (d even), `(-1)^{(d-1)/2} A_d/d!` (d odd); both cases equal
`eps_{d-1} A_d/d!` with `eps_m = (-1)^floor((m+1)/2)` (check `d mod 4`). For the second identity,
`(1-sinh z)/cosh^2 z = sech^2 z - sech z tanh z = d/dz(sech z + tanh z)`, so
`[z^m](1-sinh z)/cosh^2 z = (m+1) delta_{m+1} = (m+1) eps_m A_{m+1}/(m+1)! = eps_m A_{m+1}/m!`. ∎
[Both identities checked symbolically for `d, m <= 20`, part G.]

"""
p2 = p2.replace(old61, new61)
open('proof_part2.md','w').write(p2)
print('patched proof_part2.md')

# assemble proof_lead.md = part1 + part2 + appendix
app = r"""
## 7. Verification guide (verify_lead2.py, ALL PASS, < 15 min)

| part | statement | check | digits/range |
|---|---|---|---|
| A | D-table = binomial x zigzag, d=1..9 | exact rationals from kern_1..18 (Newton differences) | exact |
| B | u-measures m_odd, m_e, m_diag, m_oddlast; t=(1+u)/2; m_odd=1/(1-u)-u m_e | sympy radsimp == 0 | symbolic |
| C | odd-d chain: key sum = T_d - piece2 = (d!/2) e_{d+1}; even-d: key sum = -(1/2) int log^d m_e = (d!/2) e_{d+1} | quadrature dps 70 | d odd 1..9, even 2..8; 55-70 digits |
| D | closed forms T_d = d!(1-2^{-d}) zeta(d+1); piece2 = T_d - (d!/2)e_{d+1} | quadrature dps 70 | d = 1..11 odd |
| E | Master Lemma: WW = E+E-, FF = (E-E-)/s, GG = 2pi/sin - (E+E-)/s; Beta sub (4); IBP; F-series; recursion (10) | quad + nsum dps 60 | strip points s = 0.3, -0.4, 0.73, 0.2i, ...; 26-62 digits |
| F | full pipeline kernels -> c(w), w = 2..10 | quadrature dps 100 | 91-102 digits vs law |
| G | E(y) Gamma identities, duplication, log E expansion, zigzag signs, Hadamard GF form | sympy + series | symbolic / 30 digits |

All numbers printed by the scripts; no identity in this document is asserted beyond the ranges above
except those proved in writing for all d (Lemmas 2.1, 4.1, 5.1, 5.2, 5.3, Theorems 4.2, 5.4, 5.5, 6.1,
Master Lemma).
"""
full = open('proof_part1.md').read() + p2 + app
open('proof_lead.md','w').write(full)
print('wrote proof_lead.md', len(full))
with open('progress.log','a') as f:
    f.write(f"{time.strftime('%H:%M:%S')} step: assembled proof_lead.md (parts 1+2+appendix); fixed 5.5/6.1\n")
