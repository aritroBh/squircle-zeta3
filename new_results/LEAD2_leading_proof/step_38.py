import os
os.chdir('<project>/research/08_theory/LEAD2_leading_proof')
p1 = open('proof_part1.md').read()
p2 = open('proof_part2.md').read()
# check whether p2 ends with Corollary (no section 7)
print('p2 has sec7:', '## 7.' in p2, ' p1 has sec7:', '## 7.' in p1)
app = r"""
## 7. Verification guide (verify_lead2.py, ALL PASS, < 15 min)

| part | statement | check | digits/range |
|---|---|---|---|
| A | D-table = binomial x zigzag, d=1..9 | exact rationals from kern_1..18 (Newton differences) | exact |
| B | u-measures m_odd, m_e, m_diag, m_oddlast; t=(1+u)/2; m_odd=1/(1-u)-u m_e | sympy radsimp == 0 | symbolic |
| C | odd-d chain: key sum = T_d - piece2 = (d!/2) e_{d+1}; even-d: key sum = -(1/2) int log^d m_e = (d!/2) e_{d+1} | quadrature dps 70 | d odd 1..9, even 2..8; 55-70 digits |
| D | closed forms T_d = d!(1-2^{-d}) zeta(d+1); piece2 = T_d - (d!/2)e_{d+1} | quadrature dps 70 | d = 1..11 odd |
| E | Master Lemma: WW = E+E-, FF = (E-E-)/s, GG = 2pi/sin - (E+E-)/s; Beta sub (4); IBP; F-series; recursion (10) | quad + nsum dps 60 | strip points; 26-62 digits |
| F | full pipeline kernels -> c(w), w = 2..10 | quadrature dps 100 | 91-102 digits vs law |
| G | E(y) Gamma identities, duplication, log E expansion, zigzag signs, Hadamard GF form | sympy + series | symbolic / 30 digits |

All numbers printed by the scripts; no identity in this document is asserted beyond the ranges above
except those proved in writing for all d (Lemmas 2.1, 4.1, 5.1, 5.2, 5.3, Theorems 4.2, 5.4, 5.5, 6.1,
Master Lemma).
"""
full = p1 + p2 + app
open('proof_lead.md','w').write(full)
print('proof_lead.md reassembled:', len(full), 'sections:', [s for s in ['## 1.','## 2.','## 3.','## 4.','## 5.','## 6.','## 7.'] if s in full])
