# Referee report: Theorems `thm:M2` and `thm:M3` of `paper/ieee/main.tex`

Date: 2026-11-22. Independent referee. Nothing outside this folder was changed. The project scripts were
not imported or run. Every check below uses code in this folder.

## Verdicts

| Theorem | Verdict |
|---|---|
| `thm:M2` (closed form (T1) of $M_2$) | **PROVED.** The paper's sketch needs one small wording repair (R1). |
| `thm:M3` ($M_3=W_4+W_3'+T$, $W_3'$ = (T2)) | **PROVED AFTER REPAIRS.** The mathematics is right. The written proof was missing: the notes give 11 lines and the paper one sentence. The full argument is written out below (R2, R3). |

## 1. Statement check (`s1_numeric.py`, logs `s1_dps80.log`, `s1_dps160.log`)

My code takes $M_k=-\tfrac12[\varepsilon^k]\int_0^1\Phi_\varepsilon\,dt$ straight from Lemma `lem:holo`. On
$t\le\frac12$ it uses the form of $\Phi$ from the lemma's own proof, which has no cancellation, and it
expands in $\varepsilon$ with sympy. The integral is done by mpmath tanh-sinh. The closed forms (T1), (T2), $T$ and the
$W_4$ integral are copied from `main.tex` itself.

| Check (computed at 160 digits) | agreement |
|---|---|
| $M_1$ (definition) vs `thm:M1` (control) | 159.9 digits |
| $M_2$ (definition) vs (T1) | **160.2 digits** |
| $M_3$ (definition) vs $W_4$ (printed integral) + (T2) + $(\pi^2+6\lambda^2)/24$ | **159.2 digits** |
| $T$ integral vs $(\pi^2+6\lambda^2)/24=(\zeta(2)-\mathrm{Li}_2(\frac12))/2$ | 158.8 / 161 digits |
| printed decimals of $M_2$, $M_3$, $W_4$ (20 digits each) | all correct, correctly rounded |
| eq. `chi3` (side remark) | 159.5 digits (numerical only; its certificate was not part of this check) |
| $\mathrm{Li}_3(\frac12)$ formula; $(11\pm15\sqrt2)$ rewriting | exact / algebraically immediate |

Values: $M_2=-0.022161911087269399423754937739791514978\ldots$,
$M_3=0.022770521491786577456730509742327800250\ldots$, $W_4=0.795785059832240811320359819629726762547\ldots$

## 2. Proof check, step by step

Sources: `SN` = `paper/squircle_note/main.tex`, lines 166-266 (setup), 336-430 ($M_2$) and 455-466 ($M_3$). Also `main.tex`, app. `app:M2`, and `sec:method`.
Categories: (a) proved by hand in the text; (b) an exact computer-algebra identity; (c) numerical only.

| Step | Where | Cat. | My check |
|---|---|---|---|
| S0. $F$ is holomorphic; differentiate under $\int$ | `lem:holo` | (a) | Read line by line. The bounds $\lvert b_\varepsilon\rvert\le2^{1-\mathrm{Re}\,\varepsilon}\le\lvert a_\varepsilon\rvert$ and $\lvert1+\sqrt{1+z}\rvert\ge1$ give $\lvert\Phi\rvert\le2^{2-\mathrm{Re}\,\varepsilon}$. Correct. |
| S1. Euler substitution, $\log t,\log(1-t)$ split into letters | SN `lem:euler` | (a) (trivial) | Re-derived exactly (`s2_ibp.py`). |
| S2. $t$-integrands of $M_2,M_3$ (Faà di Bruno) | SN eq. M2t, `mk_integrands.py` | (b) | Re-derived exactly from $\Phi$ (`s2_ibp.py`). I did not use eq. M2t. |
| S3. Fold $v\to1/v$; integrate the poles at $\pm i$ by parts with $G(b)=0$; kernels of $W_3,T$ ($k=2$) and of $W_4,W_3',T$ ($k=3$) | SN `prop:split`, `ibp_reduce.py` | (b), with the boundary-term argument by hand (a) | **Re-derived exactly** with my own Hermite reduction. All kernels are identical to the printed ones (difference 0). The simple-pole parts at $\pm i$ are exactly 0. No weight-one or polynomial remainder is left. Every raw monomial kernel vanishes at the endpoint where its logarithms blow up, so each piece converges on its own. The boundary-term argument in the note is correct ($i\ge1$ after folding; $G(b)=0$). |
| S4. $T=(\pi^2+6\lambda^2)/24$ | SN `prop:T`; `m2_decomposition_and_tail.py` | (b), with an outline by hand | **Re-derived exactly** by a different route (my own primitives). The exact difference is 0. |
| S5. $W_3$ written exactly as 15 $\mathrm{Li}_3$ + 17 $\mathrm{Li}_2$ (split at $v=1$, polarization, primitives) | SN sec. M2.2; `m2_weight3_proof.py` | (b). The primitives are (a), checkable by differentiation. | **Re-derived independently** (`s3_weight3.py`). I get the same counts, 15 $\mathrm{Li}_3$ and 17 $\mathrm{Li}_2$, and the expression matches quadrature to $3\cdot10^{-58}$ at 60 digits. I checked the primitives by differentiation in all three regions ($10^{-51}$) and their limit $-2\zeta(3)$ at $s\to1^\pm$. |
| S6. Functional equations (I3), (T3), (D3), (R), (I), (La), (D), (A), (K) | SN; app. `app:M2` | standard theorems, quoted | All are quoted correctly with their real domains; see §3. |
| S7. `lem:trilog` (10 lines), `lem:dilog` (6 lines) | SN; `li3_table.py` | (b) | **Re-derived exactly.** My own elimination over $\mathbb Q(\sqrt2)$ (`s4_reduce.py`) uses only instances of (I3), (T3), (D3) in the argument set, the special values, and (K) at exactly the note's four points. Each of the 16 lemma lines comes out identical, with exact difference 0. |
| S8. Assembly: $W_3+T-(\text{T1})=0$ | SN; app. `app:M2` | (b) | **Re-derived**: exact difference 0 (`s4_note.log`). |
| $M_3$: $W_3'$ = (T2) | SN sec. M3 (pointer to `m3_tails_exact.py`) | (b), text missing | **Re-derived**: 18 $\mathrm{Li}_3$ + 18 $\mathrm{Li}_2$ raw (quadrature agreement $9\cdot10^{-58}$). After reduction with the same four (K) instances, (T2) is reached with exact difference 0. |

No step of either proof is only numerical (c). The numerical agreements in the notes (45 and 60 digits)
are safeguards only, and the proof does not need them. Without (K) the reduction fails: six trilogarithm values in the orbits of $\pm b^2$ stay unreduced (`s4_noK.log`). So (K) is really needed, as the note says. Allowing every (K) and Abel instance inside the argument set gives the same result (`s4_full.log`).

## 3. Functional equations: domains and branches (`s5_funeq.log`)

* (I3) $z<0$, (T3) $0<z<1$, (D3) $\lvert w\rvert\le1$, (R) $0<z<1$, (I) $z<0$, (La) $z<1$, (D): standard, and quoted
  correctly. Each was checked at 100 random points (max residual $4\cdot10^{-50}$), and every instance actually used was checked to $10^{-80}$.
* Abel (A) on $\{x<1,y<1,x+y<1\}$: correct (300 random points, $4\cdot10^{-50}$).
* Spence-Kummer (K) in SN's real form on $D=\{0<y<1,\ x<y\}$. The extension argument in the note is valid. On $D$ the nine
  arguments are all $<1$: $x/y<1$ and $x(1-y)/(y(1-x))<1$ both follow from $x<y$; $-x(1-y)/(1-x)<1$ for $x<0$ follows from $xy<1$; and the rest are $\le0$ or clearly $<1$. The logarithms $\log y$ and $\log\frac{1-y}{1-x}$ have positive arguments. So both sides are real-analytic on the convex, hence connected, set $D$. They agree on $0<x<y<1$, so the identity theorem applies. Numerically: 300 random points of $D$ (including very negative $x$) give $2\cdot10^{-49}$, and the four instances give $\le5\cdot10^{-51}$. Outside $D$ the same real-part formula fails (e.g. $(1.5,0.5)$: residual 1.51; $(-0.5,-0.2)$: 5.88), so the domain restriction matters. The four points $(b,\frac1{\sqrt2})$, $(2-\sqrt2,\frac1{\sqrt2})$, $(\frac1{\sqrt2},2b)$ and $(-1,b^2)$ all lie in $D$. The only one off the classical range is $(-1,b^2)$, and it is covered by the extension.
* Branches: every $\mathrm{Li}_n$ argument that enters the final reduction is $<1$, so no branch cut is used. The paper's
  convention "real part on $z>1$" is never needed for these two theorems.
* One reservation: I could not consult Lewin's book. The citation "[Lewin81, Eq. (6.107)]" is checked only by the fact that the
  quoted form is numerically an identity on $D$. A reader who wants to verify it should compare against the book.

## 4. Does "same as $W_3$" carry over to $W_3'$? Yes. Here is the full argument (repair R3, text for SN sec. M3)

> **Proof of Theorem M3.** By Proposition `prop:split` ($k=3$),
> $M_3=W_4+W_3'+T$ with
> $W_3'=\int_b^a\bigl[\ell_t^2\,\Omega_{20}+\ell_t\ell_{1-t}\,\Omega_{11}\bigr]dv$,
> $\Omega_{20}=\frac{1-\sqrt2}{4v}+\frac1{4(v-a)}-\frac1{4(v+a)}-\frac1{2(v+b)}$ and
> $\Omega_{11}=\frac{1+\sqrt2}{4v}-\frac1{2(v+a)}$. Here $\ell_t=\ell_b+\ell_{-a}-\ell_0-2\lambda$ and $\ell_{1-t}=\ell_a+\ell_{-b}-\ell_0-2\lambda$,
> with $\ell_x=\log|v-x|$. Proposition `prop:T` gives $T$. Expanding gives a $\mathbb Q(\sqrt2)$-combination of
> $\int\ell_x\ell_y\,dv/(v-z)$, $\lambda\int\ell_x\,dv/(v-z)$ and $\lambda^2\int dv/(v-z)$, with $x\in\{0,\pm a,\pm b\}$ and
> $z\in\{0,a,-a,-b\}$.
> (i) *Poles off $[b,a]$* ($z=0,-a,-b$). Split $[b,a]$ at $v=1$, so that each half contains at most one letter, at an
> endpoint. Then proceed exactly as for $W_3$. For $x\ne y$, polarize $\ell_x\ell_y=\frac12(\ell_x^2+\ell_y^2-\log^2|w|)$ with
> $w=\frac{v-x}{v-y}$, putting the endpoint letter in the numerator so that $w\to0$, never $\infty$. Use
> $d\log|v-z|=\frac{dw}{w-c}-\frac{dw}{w-1}$ with $c=\frac{x-z}{y-z}$, or $-\frac{dw}{w-1}$ if $y=z$, or $\frac{dw}{w}-\frac{dw}{w-1}$ if $x=z$. Use
> $\ell_x=\log|z-x|+\log|s|$ with $s=\frac{v-x}{z-x}$ for the squares. Every term then becomes $\int(\kappa+\log|s|)^n\,ds/(s-1)$,
> $n\le2$, along a path on which $s$ stays in one of $(-\infty,0)$, $(0,1)$, $(1,\infty)$. The primitives quoted in SN apply.
> (ii) *The pole at the endpoint $v=a$* (only in $\Omega_{20}$, against $\ell_t^2$). Since
> $\ell_t(a)=\log2+\log(2a)-\log a-2\lambda=0$, we can write $\ell_t=\sum_xc_x m_x$ with $m_x=\log\bigl|\frac{v-x}{a-x}\bigr|$,
> $(c_b,c_{-a},c_0)=(1,1,-1)$. No letter $x$ equals $a$. Then $\ell_t^2/(v-a)=\sum c_xc_y\,m_xm_y/(v-a)$. For $x=y$,
> $s=\frac{v-x}{a-x}$ gives $\int\log^2|s|\,ds/(s-1)$ with $s\to1$ at $v=a$. For $x\ne y$, polarization gives
> $m_x-m_y=\log|w/c|$ with $c=w(a)=\frac{a-x}{a-y}$, and
> $\int\log^2|w/c|\bigl(\frac{dw}{w-c}-\frac{dw}{w-1}\bigr)$. In the first part $s=w/c\to1$ again. The second part is regular at $v=a$.
> At $s\to1$ (from either side) the primitive $\log^2s\log(1-s)+2\log s\,\mathrm{Li}_2(s)-2\mathrm{Li}_3(s)$ (or its
> $s>1$ form) tends to $-2\zeta(3)$. All coefficients $\kappa$ vanish there, so no $\log|1-s|$ divergence occurs. Every term is
> integrable at $v=a$, since $m_x$ vanishes linearly there.
> (iii) All arguments are cross-ratios of $\{0,\pm a,\pm b,1,\infty\}$, hence $S$-units of $\mathbb Z[\sqrt2,\frac12]$. This writes
> $W_3'$ exactly with 18 values of $\mathrm{Li}_3$ and 18 of $\mathrm{Li}_2$, at the arguments that occur for $W_3$, together with $\mathrm{Li}_3$ at $-a^2$, $-\sqrt2$, $b^2$ and $\mathrm{Li}_2(-\sqrt2)$.
> (iv) The $\mathrm{Li}_2$ values reduce to $\mathrm{Li}_2(b)$ by Lemma `lem:dilog`, (R), (La) and (I). The $\mathrm{Li}_3$ values reduce to
> $\mathrm{Li}_3(\pm b)$, $\mathrm{Li}_3(1/\sqrt2)$ and $\zeta(3)$ by Lemma `lem:trilog`, together with (I3) at $-a^2$ and at $-\sqrt2$ (which send them to $-b^2$ and $-1/\sqrt2$) and (D3) at $b$ (for $b^2$). No
> Kummer instance beyond the four used for $W_3$ is needed. The result is exactly (T2). $\square$

I carried out exactly this argument in `s3_weight3.py` and `s4_reduce.py`. The raw expression agrees with quadrature to
$9\cdot10^{-58}$, and the exact difference from (T2) is $0$.

## 5. Repairs and sentences to change

**R1 (main.tex, app. `app:M2`, weight-two part).** The paper says "With Lemma~\ref{lem:dilog} all the $\Li_2(b)$ terms cancel". The
list of dilogarithms just before it includes $\mathrm{Li}_2(2b)$, $\mathrm{Li}_2(-a^2)$ and $\mathrm{Li}_2(1/\sqrt2)$, and Lemma `lem:dilog` does not
cover them. Suggested text: "With Lemma~\ref{lem:dilog}, reflection at $2b$
($\Li_2(2b)=\frac{\pi^2}6+2\lambda L-2L^2-\Li_2(b^2)$), Landen at $-a$ and inversion at $-a$ and $-a^2$, all the $\Li_2(b)$ terms cancel
and $T=(\pi^2+6\lambda^2)/24$."

**R2 (main.tex, app. `app:M2`, last sentence).** Replace "The same procedure gives $W_3'$; the one new feature is a pole at the
endpoint $v=a$, where the primitive takes the value $-2\zeta(3)$." with:
"The same procedure gives $W_3'=\int_b^a[\ell_1^2\Omega_{20}+\ell_1\ell_2\Omega_{11}]dv$, with
$\Omega_{20}=\frac{1-\sqrt2}{4v}+\frac1{4(v-a)}-\frac1{4(v+a)}-\frac1{2(v+b)}$ and $\Omega_{11}=\frac{1+\sqrt2}{4v}-\frac1{2(v+a)}$. The one new
feature is the pole of $\Omega_{20}$ at the endpoint $v=a$. Since $\ell_1(a)=0$ we write $\ell_1=\sum_xc_x\log|\frac{v-x}{a-x}|$. After
polarization the singular terms become $\int\log^2|s|\,ds/(s-1)$ with $s\to1$, where the primitive takes the value $-2\zeta(3)$. The
same four Kummer instances suffice." Without the integrand for $W_3'$, the paper's statement "$M_3=W_4+W_3'+T$" defines $W_3'$ only
through its closed form, and the sketch cannot be followed.

**R3 (SN sec. M3).** Replace the 11-line proof with the argument in §4.

**Other sentences:**
* main.tex, after `thm:M2`: "a 290-digit search finds no way to write $M_2$ with $\zeta(3)$ and products of lower-weight constants
  alone, so trilogarithms are really needed (evidence, not proof)." The words "really needed" claim more than the parenthesis allows.
  Suggest: "so trilogarithms appear to be needed (evidence, not proof)."
* `DAGGER_INVENTORY.md`, row 2, says the paper gives "statement only" for `thm:M3`. In fact app. `app:M2` is titled "Proofs of Theorems M2 and
  M3" and has one sentence on $W_3'$. That sentence is too thin (R2).
* Statements of `thm:M2` and `thm:M3` as worded: accurate. "A combination of weight three of $\zeta(3)$, $\mathrm{Li}_3(\pm b)$,
  $\mathrm{Li}_3(1/\sqrt2)$ and products of lower-weight constants" is correct: the only other term is $\mathrm{Li}_2(b)$ times a logarithm, which is
  a product. All printed decimals are correct.
* The app. `app:M2` claims "15 trilogarithms and 17 dilogarithms", "in the real form valid on $\{0<y<1,\ x<y\}$", and "the difference … is
  exactly zero". All three were confirmed independently.

## Files
`s1_numeric.py` (statement check), `s2_ibp.py` (exact weight split, k=1,2,3), `q2.py` (exact $\mathbb Q(\sqrt2)$ and S-unit logs),
`s3_weight3.py` (own exact evaluation of $T$, $W_3$, $W_3'$), `s4_reduce.py` (own exact reduction; modes `note`, `full`, `noK`),
`s5_funeq_checks.py` (primitives and functional equations), and the logs `*.log`.
Run with `<project>/.venv/bin/python`; each script takes under 25 s.
