# W2SPLIT: the reduction splits every M_k (2026-11-23)

**Status: PROVED for every k >= 2, by hand.** No simple pole at +-i survives at any degree, at any level.
Exact check of every step of the lemma: k <= 40 (all 440 keys at k = 40). This is a new argument and needs a referee.

## The precise lemma

Notation of `paper/ieee/proofs/sec_g1_wt2.tex` (Definition g1:def:red): kernels W^(k)_rs, primitives G^(k)_rs = J(W),
H(W) = principal parts of order >= 2 at +-i. Put iota(v) = -1/v and I(W)(v) = W(-1/v)/v^2 (pull-back of W dv).

**Lemma (parity).** For every k >= 2 and every key (r,s): W^(k)_rs + I(W^(k)_rs) has no pole at +-i; for d >= 2,
G^(k)_rs(v) + G^(k)_rs(-1/v) is constant.

**Consequence (what the caveat asks).** Res_{v=+-i} W^(k)_rs = 0 for every key, so the simple kernel
S_rs = W_rs - H(W_rs) = W_rs - G_rs' has no pole at +-i. Moreover (Lemma "real letters")
S_rs = sum over x in {0, a, -a, -b} of c_x/(v-x), c_x in Q(sqrt2), no pole at b, and (Lemma "splits")
M_k = sum over keys with 1 <= r+s <= k of int_b^a S_rs l1^r l2^s dv, every integral absolutely convergent,
the d = 1 term being int_b^a R_w2 l1 dv = T.

## Idea of the proof (three lines)

1. v -> t(v) is two-to-one; the other preimage of t(v) is iota(v) = -1/v, which fixes +-i. t is iota-invariant and
   q = (v^2+1)/(2 sqrt2 v) is iota-odd, so the square-root part of the integrand is iota-odd and the rest (E_rs, multiples
   of D_1, D_2) is iota-even and regular at +-i.
2. The reduction respects this sign: D_j dv = d l_j is invariant, I commutes with Fold and with H, so by descending
   induction every W is iota-odd modulo functions regular at +-i, and every G is iota-odd up to a constant.
3. The residue at a fixed point is invariant under I (I(1/(v-i)) = 1/(v-i) - 1/v), so an iota-odd W has
   2 Res_i W = 0. Hermite reduction then leaves nothing at +-i.

The generating-function view (remark in the fragment): in the t-plane the radicand t^(2e-2)+(1-t)^(2e-2) has, for small
e != 0, one simple branch point near (1+i)/2; the Euler substitution is the double cover that uniformises it, and a
1-form that is odd under the sheet exchange has no residue at the ramification point. This is why the statement holds
for all levels at once.

What we use from the supplement: Lemma g1:lem:init (closed form of the initial data; we only use its shape
P_rs + E_rs, not the constants), Lemma g1:lem:A (pole bound), Lemma g1:lem:kernel (W_10 = R_w2), Lemma g1:lem:euler,
the folding rule, and Lemma g1:lem:T. Nothing else.

## What is by hand and what is computation

By hand (in `paper/ieee/proofs/sec_g1_wt2_split.tex`): Lemmas g1:lem:inv (properties of I), g1:lem:par (parity),
g1:lem:res (no residue), g1:lem:real (kernels are simple at 0, a, -a, -b; no pole at b; this includes the cancellation
of the poles of P_{d0} and E_{d0} at b), g1:lem:split (the identity, convergence and vanishing boundary terms).
No step depends on a script.

Exact computation (checks, not proofs), `split_check.py`, own code (python-flint 0.9, exact arithmetic in Q(sqrt2)(v)),
independent of `paper/ieee/proofs/scripts/` and of the notes' reducers. Initial data are built from the two-variable
generating function Psi(v;X,Y) = -1/2 (t^-1 e^X + (1-t)^-1 e^Y - S sqrt(1+u)) dt/dv with Stirling numbers
([X^r](e^{2X}-1)^a = 2^r a! S(r,a)/r!), not from the formula of Lemma g1:lem:init. For k = 2..40 and every key:
- (P) W + I(W) has no pole at +-i;
- (R) residue of W at +-i is 0 (Hermite reduction at v^2+1 leaves residual numerator 0);
- (L) S = W - H(W) has denominator dividing v (v^2+2v-1)(v^2-2v-1), is O(1/v) at infinity, no pole at b;
- (G) G + G(-1/v) is constant (d >= 2);
- (K) W_10 = R_w2.
All True for k <= 40 (`split_check_k40.log`, 33 s). Controls in the same log show the checks can fail:
1/(v^2+1) and 1/(v^2+1)^2 leave a simple part; with Fold in place of I the parity fails for 7 of 8 level-4 keys;
the level-4 kernels have (v^2+1)^6, (v^2+1)^4, (v^2+1)^2 in their denominators (so the reduction is not vacuous).

Numerical safeguard (`numeric_check.py`, `numeric_check.log`, 30 digits): for k = 2..5 the kernels agree with a direct
Taylor expansion of Phi_e in e (difference ~1e-32), and the sum of the split pieces int S_rs l1^r l2^s equals
M_k = -1/2 int_0^1 [e^k] Phi_e dt computed directly (differences 3e-27 .. 7e-24, limited by the quadrature).
M_2 = -0.0221619110872693994..., as in Theorem 2.

Rerun: `cd new_results/W2SPLIT_reduction_splits; ../../../.venv/bin/python -B split_check.py 40 controls`
and `../../../.venv/bin/python -B numeric_check.py 5` (about 4 min).

## The LaTeX fragment

`paper/ieee/proofs/sec_g1_wt2_split.tex`: one subsection (label g1:sub:wt2:split), lemmas g1:lem:inv, g1:lem:par,
g1:lem:res, g1:lem:real, g1:lem:split and a remark. To be inserted in `sec_g1_wt2.tex` after Section g1:sub:wt2:proof
(the proof of Theorem 6) and before "Machine steps". Uses only macros of `macros_g1.tex` and \Res of the wrapper;
\mathcal I is written inline (no new macro). Test build (copy of M2, M3, wt2 + fragment, tectonic -C --keep-logs):
0 undefined references, no overfull box above 2pt from the fragment.

## OLD/NEW text for sec_g1_wt2.tex (not applied)

1. Line 14.
   OLD: `For $k=2,3$ Sections~\ref{g1:sec:M2} and~\ref{g1:sec:M3} show that the reduction really splits $M_k$ into pieces, with the weight-two piece $T$. For $k\ge4$ the reduction is used here as the paper uses it; the statement that it splits $M_k$ (no simple pole at $\pm i$ survives at degrees $\ge2$) is not part of Theorem~6 and is not proved here (the paper reports an exact check for $k\le12$).`
   NEW: `For $k=2,3$ Sections~\ref{g1:sec:M2} and~\ref{g1:sec:M3} show that the reduction really splits $M_k$ into pieces, with the weight-two piece $T$. Section~\ref{g1:sub:wt2:split} proves this for every $k\ge2$: no simple pole at $\pm i$ survives at any degree, and $M_k$ is a finite sum of convergent integrals $\int_b^a\ell_1^r\ell_2^s\,\frac{dv}{v-x}$, $x\in\{0,a,-a,-b\}$, whose degree-one term is $T$ (Lemma~\ref{g1:lem:split}). Theorem~6 itself does not need this.`

2. Line 148 (Machine steps).
   OLD: `\emph{What is by hand:} Lemmas~\ref{g1:lem:init}, \ref{g1:lem:A}, \ref{g1:lem:B}, \ref{g1:lem:I}, \ref{g1:lem:D}, \ref{g1:lem:kernel} and the theorem;`
   NEW: `\emph{What is by hand:} Lemmas~\ref{g1:lem:init}, \ref{g1:lem:A}, \ref{g1:lem:B}, \ref{g1:lem:I}, \ref{g1:lem:D}, \ref{g1:lem:kernel} and the theorem, and Lemmas~\ref{g1:lem:inv} to~\ref{g1:lem:split} on the splitting;`

3. Machine-steps table, new row after the row of `g1\_wt2\_check.py 7`:
   NEW: `\texttt{split\_check.py 40 controls} (notes: \path{W2SPLIT_reduction_splits/}) & the reduction of Definition~\ref{g1:def:red} for $k=2,\dots,40$ (python-flint, exact over $\Qs$; initial data from the two-variable generating function by Stirling numbers, independent of \texttt{g1\_wt2\_check.py}): for every key, $W+\mathcal I(W)$ has no pole at $\pm i$, $\Res_{\pm i}W=0$, $G(v)+G(-1/v)$ is constant, the simple kernel has poles only at $0,\pm a,-b$, and $W^{(k)}_{10}=\Rwt$; negative controls. Every line prints \texttt{True}. & 33 s\\`
   (If the planner prefers the scripts folder: copy `split_check.py` to `paper/ieee/proofs/scripts/g1_wt2_split_check.py` and change the name in this row.)

4. Line 175 (Checks, "Exact, finite").
   OLD: `and finds $W^{(k)}_{10}=\Rwt$, the leading coefficients of Lemma~\ref{g1:lem:D} for every key, and no residual simple pole at $\pm i$ at any degree. The notes report the same for $k\le15$ (weight two) and, for the recursion, $k\le12$; the paper's remark is the same. These runs are checks, not proofs.`
   NEW: `and finds $W^{(k)}_{10}=\Rwt$, the leading coefficients of Lemma~\ref{g1:lem:D} for every key, and no residual simple pole at $\pm i$ at any degree. The notes report the same for $k\le15$ (weight two) and, for the recursion, $k\le12$; the paper's remark is the same. These runs are checks, not proofs; the absence of a simple pole is proved for every $k$ in Lemma~\ref{g1:lem:res}. An independent script (notes: \path{W2SPLIT_reduction_splits/split_check.py}) checks Lemmas~\ref{g1:lem:par}, \ref{g1:lem:res}, \ref{g1:lem:real} at every key for $k\le40$. \emph{Numerical:} for $k=2,\dots,5$ the sum of the integrals of Lemma~\ref{g1:lem:split} agrees with $M_k=-\frac12\int_0^1\varphi_k\,dt$, computed directly, to $23$ digits or more (30-digit quadrature).`

5. Line 177 ("Which report confirmed the proof"), append:
   NEW sentence: `Section~\ref{g1:sub:wt2:split} is new and has not yet been re-checked by an independent pass.`

## Consequence for main.tex (for the planner; not applied)

main.tex, the weights 4..7 paragraph ("polynomial in $k$ of degree"), says for w = 4, 5, 6, 7 that "Exact checks at
k = w, ... show that their simple-pole coefficients at these points vanish". Lemma g1:lem:res proves that these
simple-pole coefficients at +-i vanish for every k and every degree, so that part of each claim no longer rests on the
checks. The other checks in that paragraph (that k = w lies on the same polynomial, and the tower law itself) are not
touched by this result.

## Files

- `split_check.py`, `split_check_k40.log` (exact, k <= 40, with controls)
- `numeric_check.py`, `numeric_check.log` (numerical safeguard, k <= 5)
- `result.md` (this file)
- `paper/ieee/proofs/sec_g1_wt2_split.tex` (LaTeX fragment, new file)
