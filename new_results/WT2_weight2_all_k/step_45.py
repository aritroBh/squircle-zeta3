result = r'''# WT2: Conjecture 3 — the weight-two part of M_k is (pi^2+6 lambda^2)/24 for all k >= 2

## Status: exact checks to k = 15 (kernel identity), structure identified, proof in progress

## 1. Reproduction and extension of the exact check (Task 1)

The reduced M_k is a sum over log-degrees; the weight-two part is the degree-1 piece
`int_b^a K_k(v) l1 dv` with l1 = log t. The conjecture is equivalent to the **kernel identity**

    K_k(v) = Rw2(v) = 1/(2v) + 1/(2(v-a)) - 1/(v+a) - 1/(2(v+b))     for all k >= 2,

since T = int_b^a Rw2(v) log t(v) dv = (pi^2+6 lambda^2)/24 is proved in the paper (Appendix, dilog identities).

Exact results (kernel - Rw2 == 0 in exact arithmetic):

| k | method | time (build+reduce) | kernel == Rw2 |
|---|--------|---------------------|----------------|
| 2..5 | original `ibp_reduce.reduce_k` (sympy apart/integrate over Q(sqrt2,i)) | 2s, 6.6s, 22s, 73s | YES (reproduced) |
| 2,3 | `fastred.reduce_k_fast`, cross-checked key-by-key vs original on ALL keys | — | YES (identical) |
| 4..8 | fast reducer on converted `v_integrand` | 0.1s, 0.3s, 0.7s, 1.3s, 2.7s (reduce) | YES |
| 2..15 | closed-form integrand `build_integrand_pf` + fast reducer | see table | YES |

No failures. The identity `K_k = Rw2` holds exactly for every 2 <= k <= 15.

### Timings for the full exact check (closed-form integrand + exact PF reducer)

| k | build (s) | reduce (s) |
|---|-----------|------------|
| 9  | 5.5  | 4.1  |
| 10 | 7.6  | 6.5  |
| 11 | 11.0 | 9.4  |
| 12 | 15.9 | 15.7 |
| 13 | 30.0 | 19.3 |
| 14 | 33.7 | 30.3 |
| 15 | 35.5 | 56.6 |

(higher k running)

## 2. The machinery (all exact, no sympy in the certified path)

* `fastred.py`: exact arithmetic in Q(sqrt2,i) as 4 Fractions; partial-fraction
  representation of all integrands (simple real poles in {0,+-a,+-b}, principal parts
  at +-i, no polynomial part); fold v->1/v, split simple/higher, rational antiderivative
  of the higher part, evaluation at b, multiplication by DL1/DL2 -- all in closed form.
  `pf_from_poly` verifies the denominator shape `lead*(v^2+1)^m * prod (v-r)`.
* `genintegrand.py`: closed form of the M_k integrand. The coefficient of l1^m l2^{k-m} in
  the M_k integrand is
      c_m(t) = q S_m(t) / (2 t (1-t) q2^k)  -  delta_{m,k}/(2 k! t)  -  delta_{m,0}/(2 k! (1-t)),
  q2 = t^2+(1-t)^2, q = sqrt(q2), and
      S_m(t) = sum_{j=1}^k binom(1/2,j) sum_{s=0}^j binom(j,s) w(s,m) w(j-s,k-m)
               t^{2(j-s)} (1-t)^{2s} q2^{k-j},     w(s,m) = s! S(m,s) 2^m/m!  (S = Stirling 2nd).
  Verified: matches `mk_integrands.t_integrand` for k = 2..5 (after Q = sqrt(q2)) and
  `v_integrand` for k = 2..8 (after conversion to v), in exact arithmetic.
* `buildintegrand.py`: conversion to v via t=(v^2+2v-1)/(4v) (exact polynomial arithmetic),
  folding at the t level, A/B parts.

## 3. Structure / proof (Task 2) — in progress
(to be filled)

## 4. Paragraph for the paper (Task 4) — in progress
'''
with open('result.md','w') as f:
    f.write(result)
print("result.md written")
