
## 4. Even d (odd weights): the LEAD1 argument, refereed

**Lemma 4.1 (parity collapse).** For even `d`,  `K_d = -(1/2) int_{-1}^1 (L1-L2)^d m_e(u) du`.

**Proof.** Write `T_j := (-1)^{j+1} binom(d,j) M(d-j,j; m_e)`. Since `m_e` is even, under `u -> -u`
(`L1 <-> L2`) `M(d-j,j;m_e) = M(j,d-j;m_e)`; and since `d` is even, `(-1)^{d-j+1} binom(d,d-j) =
(-1)^{j+1} binom(d,j)`. Hence `T_j = T_{d-j}` and `sum_{j<d/2} T_j = (sum_{j=0}^d T_j - T_{d/2})/2`.
The full sum is `sum_j (-1)^{j+1} binom(d,j) M(d-j,j;m_e) = -int sum_j binom(d,j) L1^{d-j}(-L2)^j m_e du
= -int (L1-L2)^d m_e du`. In (3) the middle key carries `m_diag`, not `m_e`, so

    K_d = sum_{j<d/2} T_j + (-1)^{d/2+1} binom(d,d/2) M(d/2,d/2; m_diag)
        = -(1/2) int (L1-L2)^d m_e du + (-1)^{d/2+1} binom(d,d/2) M(d/2,d/2; m_diag - m_e/2).

Since `m_diag - m_e/2 = -u m_e/2` is odd in `u` while `L1^{d/2} L2^{d/2}` is even, the correction
integral vanishes.  ∎  [Collapse checked to 60-70 digits for d = 2,4,6,8, part C.]

**Theorem 4.2.** For even `d`, `K_d = (d!/2) e_{d+1}`.

**Proof.** Lemma 4.1 gives `K_d = -(1/2) int (L1-L2)^d m_e du
= -(d!/2) [z^d] int_{-1}^1 e^{z(L1-L2)} m_e(u) du`. With `u = (1-w)/(1+w)` (so `du = -2/(1+w)^2 dw`,
`1-u^2 = 4w/(1+w)^2`, `sqrt2/S = (1+w)/sqrt(1+w^2)`) one gets, for any `f`,

    int_{-1}^1 f(u) m_e(u) du = -(1/2) int_0^inf f((1-w)/(1+w)) h(w) dw,                    (5)

hence `int_{-1}^1 ((1+u)/(1-u))^z m_e du = -(1/2) FF(z)` on the strip. Therefore for even `d`

    K_d/d! = -(1/2) [z^d] (-(1/2) FF(z)) = (1/4) [z^d] FF(z)
           = (1/4) [z^d] (E(z)-E(-z))/z        [Master Lemma (ii)]
           = (1/4) * 2 e_{d+1} = e_{d+1}/2,

because `(E(z)-E(-z))/z = 2 sum_{m odd} e_m z^{m-1}` and `m = d+1` is odd.  ∎

**Referee notes on the LEAD1 version (genfun.md 3.3).** Re-derived independently:
(a) the collapse (Lemma 4.1) — OK;
(b) the `w`-substitution and the inversion symmetry `h(1/w) = w^2 h(w)` reducing `int_0^inf` to
`2 int_0^1` — OK (subsumed in (5) + Step 3);
(c) `F(s) = int_0^1 w^{-s} h dw` converges for `Re s < 1`, its partial-fraction series converges
absolutely off the poles, termwise integration justified by uniform convergence of the alternating
binomial series — OK, with the printed erratum (second sum starts at `n = 0`), Remark in Section 3;
(d) `(F(s)+F(-s))/2 = (E(s)-E(-s))/(2s)` for `-1 < Re s < 1` — OK (this is Master Lemma (ii);
the proof via `E = (1/2)B`, split at `t = 1/2`, substitutions `t = w^2/(1+w^2)`, and IBP with
boundary vanishing exactly on `-1 < Re s < 1` re-checked line by line);
(e) coefficient extraction `[y^d]` at `y = 0` — OK (`F` analytic at 0, strip contains 0).
Verdict: **OK** (one documentation typo in the printed series; no mathematical gap).

## 5. Odd d (even weights): the new proof

Throughout, `d` is odd, `mid = (d-1)/2`.

### 5.1 Decomposition of the key sum

By (2), `m_odd = 1/(1-u) - u m_e`, and `m_oddlast = m_odd - 1/(2(1-u))`. Hence (3) splits as

    K_d = [ sum_{j<=mid} (-1)^{j+1} binom(d,j) M(d-j,j; 1/(1-u))
            - (1/2)(-1)^{mid+1} binom(d,mid) M(mid+1,mid; 1/(1-u)) ]                       (6a)
          - sum_{j<=mid} (-1)^{j+1} binom(d,j) M(d-j,j; u m_e).                             (6b)

Call the bracket `T_d` and the last line `-(u m_e part)`.

**Lemma 5.1 (the `u m_e` part collapses; PROVED).**
`sum_{j<=mid} (-1)^{j+1} binom(d,j) M(d-j,j; u m_e) = -(1/2) int_{-1}^1 (L1-L2)^d u m_e(u) du`,
and its (odd-`d`) generating function is `sum_{d odd} (...) z^d/d! = GG(z)/4`. Explicitly,

    (u m_e part)/d! = (1/4)[z^d] GG(z) = (1-2^{-d}) zeta(d+1) - e_{d+1}/2.                (7)

**Proof.** `u m_e(u)` is odd, so `M(d-j,j; u m_e) = -M(j,d-j; u m_e)`; and for `d` odd the
coefficient satisfies `(-1)^{d-j+1} binom(d,d-j) = (-1)^j binom(d,j)`. The two sign flips cancel:
`b_j := (-1)^{j+1} binom(d,j) M(d-j,j;u m_e)` obeys `b_{d-j} = b_j`. Since `d` is odd there is no
fixed point, so `sum_{j<=mid} b_j = (1/2) sum_{j=0}^d b_j
= -(1/2) int sum_j binom(d,j) L1^{d-j}(-L2)^j u m_e du = -(1/2) int (L1-L2)^d u m_e du`. Applying (5) with the extra factor
`u = (1-w)/(1+w)`: `int_{-1}^1 f(u) u m_e du = -(1/2) int_0^inf f((1-w)/(1+w)) g1(w) dw`, hence
`int_{-1}^1 ((1+u)/(1-u))^z u m_e du = -(1/2) GG(z)`. Extracting the odd powers of `z`
(`-(1/2) int sinh(z(L1-L2)) u m_e du = -(1/4)[-(1/2)GG(z) + (1/2)GG(-z)] = GG(z)/4`, using that
`GG` is odd by Master Lemma (iii)) gives the GF `GG(z)/4`. The explicit form is
`[z^d] GG(z) = 2 [z^d](pi/sin pi z) - [z^d](E(z)+E(-z))/z = 4(1-2^{-d}) zeta(d+1) - 2 e_{d+1}`
(odd `d`), using `[z^d] pi/sin(pi z) = 2(1-2^{-d}) zeta(d+1)` (from
`pi/sin(pi z) = 1/z + sum_{n>=1} (-1)^{n+1} (2^{2n}-2) B_{2n} pi^{2n}/(2n)! z^{2n-1}` and
`zeta(2n) = (2 pi)^{2n} |B_{2n}|/(2 (2n)!)`). ∎  [(7) verified to 51-72 digits for d = 1..11, part D.]

### 5.2 The Bose-Einstein-log bracket T_d

Define the convergent integrals (`p >= 1`)

    A_{p,q} := M(p,q; 1/(1-u)) = int_{-1}^1 L1^p L2^q/(1-u) du
             = (-1)^p int_0^inf r^p log^q(1-e^{-r})/(e^r-1) dr                            (8)

(the second form via `(1+u)/2 = e^{-r}`). Note `M(mid+1,mid; 1/(2(1-u))) = A_{mid+1,mid}/2`, so

    T_d = sum_{j=0}^{mid} (-1)^{j+1} binom(d,j) A_{d-j,j} - (1/2)(-1)^{mid+1} binom(d,mid) A_{mid+1,mid}.   (9)

**Lemma 5.2 (two-parameter generating function; PROVED).** For `Re y > 0`, `Re x > -1`,

    H(x,y) := sum_{p,q} (A_{p,q}/(p! q!)) x^p y^q = int_{-1}^1 e^{x L1 + y L2}/(1-u) du = B(1+x, y),

where the coefficient identity `A_{p,q}/(p!q!) = [x^p y^q] B(1+x,y)` holds rigorously for all `p >= 1`,
`q >= 0` (Laurent coefficient in `y`; `B(1+x,y)` has a simple pole `1/y`).

**Proof.** The integral equals `B(1+x,y)` by `u = 2s-1` (Section 3, Step 1 pattern) or by
`s = e^{-r}` in (8): `int_0^inf e^{-xr}(1-e^{-r})^y/(e^r-1) dr = int_0^1 s^x (1-s)^{y-1} ds = B(1+x,y)`.
For `p >= 1`, `a_p(y) := (-1)^p int_0^inf r^p (1-e^{-r})^{y-1} e^{-r} dr` is analytic on `Re y > -1`
(integrand `O(r^{p+Re y-1})` at 0, exponentially small at infinity) and `A_{p,q} = a_p^{(q)}(0)`;
while `[x^p] B(1+x,y) = (1/p!) int_0^1 (log s)^p (1-s)^{y-1} ds = a_p(y)/p!` by `s = e^{-r}`. So the
Taylor coefficients in `x` agree as analytic functions of `y` near `0`, and the `y`-Laurent
coefficients agree.  ∎

**Lemma 5.3 (functional equation and recursion; PROVED).** `y B(1+x,y) = x B(1+y,x)` (both sides
equal `P(x,y) := Gamma(1+x)Gamma(1+y)/Gamma(1+x+y)`, symmetric, analytic at 0, `P(0,0) = 1`).
Hence with `aa_{p,q} := A_{p,q}/(p!q!)`,

    aa_{p,q} = aa_{q+1, p-1}     (p >= 1, q >= 0).                                        (10)

**Proof.** `aa_{p,q} = [x^p y^q] B(1+x,y) = [x^p y^{q+1}] P`; `aa_{q+1,p-1} = [x^{q+1} y^p] P`;
equal by symmetry of `P`. ∎  [Recursion checked to 58 digits for p = 1..5, q = 0..4, part E.]

**Theorem 5.4 (the A-bracket; PROVED).** For odd `d`,

    T_d/d!  =  (1/2) [z^d] (pi/sin(pi z))  =  (1 - 2^{-d}) zeta(d+1).                     (11)

**Proof.** Work on the finite set `S_d = {(p,q): p >= 1, q >= 0, p+q = d}` with
`f(p,q) := (-1)^{q+1} aa_{p,q}`. The map `sigma(p,q) = (q+1, p-1)` is a bijection of `S_d`
(it is an involution: `sigma^2 = id`), preserving `f`: `aa` is preserved by (10), and
`(-1)^{q'+1} = (-1)^{p-1+1} = (-1)^p = (-1)^{q+1}` because `p+q = d` is odd. Its unique fixed point
is `D_* = ((d+1)/2, (d-1)/2)`. With `A_* = {(p,q) in S_d: p > q}` one checks
`sigma(A_*) = (S_d \ A_*) ∪ D_*` (cardinality check: `|S_d| = d`, `|A_*| = (d+1)/2`, and
`sigma(p,q) = D_*` iff `(p,q) = D_*`). Hence `2 sum_{A_*} f = sum_{S_d} f + f(D_*)`, and

    T_d/d! = sum_{A_*} f - (1/2) f(D_*) = (1/2) sum_{S_d} f.

Since `aa_{0,d} = [y^d] B(1,y) = [y^d] (1/y) = 0` for `d >= 1`, the sum extends to `p >= 0`:
`sum_{S_d} f = sum_{p+q = d} (-1)^{q+1} aa_{p,q} = -[z^d] B(1+z, -z)`. Finally

    B(1+z, -z) = Gamma(1+z) Gamma(-z) / Gamma(1) = z Gamma(z) * Gamma(1-z)/(-z) = -pi/sin(pi z)

(reflection formula), so `T_d/d! = (1/2)[z^d](pi/sin pi z)`. The Bernoulli expansion
`pi/sin(pi z) = 1/z + sum_{n>=1} (-1)^{n+1} (2^{2n}-2) B_{2n} pi^{2n}/(2n)! z^{2n-1}` with
`pi^{2n} |B_{2n}|/(2n)! = 2^{1-2n} zeta(2n)` gives, for `d = 2n-1`,
`[z^d] pi/sin(pi z) = 2(1-2^{-d}) zeta(d+1)`. ∎  [(11) verified to 63-72 digits for d = 1..11, part D.]

### 5.3 Assembly for odd d

**Theorem 5.5.** For odd `d`, `K_d = (d!/2) e_{d+1}`.

**Proof.** By (6), `K_d = T_d - (u m_e part)`. By Theorem 5.4, `T_d/d! = (1/2)[z^d](pi/sin pi z)`.
By Lemma 5.1 and Master Lemma (iii),

    (u m_e part)/d! = (1/4)[z^d] GG(z) = (1/4)( 2[z^d](pi/sin pi z) - [z^{d+1}] WW(z) ).

The `pi/sin(pi z)` contributions cancel exactly:

    K_d/d! = T_d/d! - (u m_e part)/d!
           = (1/2)[z^d](pi/sin pi z) - (1/2)[z^d](pi/sin pi z) + (1/4)[z^{d+1}] WW(z)
           = (1/4)[z^{d+1}](E(z)+E(-z))
           = (1/4)(e_{d+1} + (-1)^{d+1} e_{d+1}) = e_{d+1}/2,

since `d+1` is even. ∎  [Full chain key sum = T_d - piece2 = (d!/2) e_{d+1} verified to 51-72 digits
for d = 1..11, part C; T and piece2 closed forms to d = 11, part D.]

**Remark.** The finite parts that blocked LEAD1 are the identical `(1-2^{-d}) zeta(d+1) d!` terms in
`T_d` and in the `u m_e` part; the involution (10) exhibits the first as a single diagonal Beta value,
and the strip representation of `GG` exhibits the second as the `2 pi/sin(pi z)` pole term. They are
literally the same obstruction seen from the two sides, and they cancel.

## 6. Assembly: the law and its generating-function forms

**Lemma 6.1 (zigzag signs).** `delta_d := [x^d](sech x + tanh x) = eps_{d-1} A_d/d!` for all `d`,
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

**Theorem 6.2 (the law).** Given the D-table of Section 1 for all `d`, for every `w >= 2`:

    c(w) = eps_{w-2} A_{w-1}/(2 (w-2)!) * e_w,   e_w = [y^w] E(y).

**Proof.** `c(d+1) = delta_d K_d/(d-1)!` (definition of `K_d`) `= delta_d (d!/2) e_{d+1}/(d-1)!`
(Theorems 4.2 and 5.5) `= delta_d (d/2) e_{d+1}`. With `delta_d = eps_{d-1} A_d/d!` (Lemma 6.1) the
prefactor is `eps_{d-1} A_d/(2 (d-1)!)`, and `w = d+1` gives the statement. ∎

**Corollary 6.3 (generating functions).**
  (a) `F(u) := sum_{w>=2} c(w) u^w = (1/2)[ u^2 (1-sinh u)/cosh^2 u  ⊙  E(u) ]`
      (Hadamard = coefficientwise product), since `eps_{w-2} A_{w-1}/(w-2)! =
      [u^{w-2}] ((1-sinh u)/cosh^2 u) = [u^w] (u^2 (1-sinh u)/cosh^2 u)` by Lemma 6.1.
  (b) The reduced GF is a pure Gamma product:
      `sum_{w>=2} [2(w-2)!/(eps_{w-2} A_{w-1})] c(w) u^w = E(u) - 1 - u log 2`.
  (c) `E(y) = 2^y Gamma(1-y) Gamma(1+y/2)/Gamma(1-y/2)`: Legendre duplication
      `Gamma(z)Gamma(z+1/2) = 2^{1-2z} sqrt(pi) Gamma(2z)` at `z = 1/2 - y/2` gives
      `Gamma(1/2-y/2) = 2^y sqrt(pi) Gamma(1-y)/Gamma(1-y/2)`, substitute in the definition of `E`.
      Hence with `log Gamma(1+t) = -gamma t + sum_{n>=2} (-1)^n zeta(n) t^n/n`,
      `log E = y log 2 + sum_{n>=2} s_n y^n/n`, `s_n = zeta(n)` (n even), `s_n = (1-2^{1-n}) zeta(n)`
      (n odd).  [All checked symbolically/numerically, part G.]
