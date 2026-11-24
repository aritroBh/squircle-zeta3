# REPORT_G8B: Proposition 5 at the level of functions

Files (all new): `sec_g8_M3func.tex` (281 lines), `macros_g8b.tex` (\providecommand only), `wrap_g8b.tex`,
`scripts/g8b_check.py`, `scripts/g8b_check.out`. No existing file was edited.

## Theorem as stated (sec_g8_M3func.tex:26, "Theorem G8B")
For 0<alpha<pi, c = cos(alpha/2), M_3 as a function of c in (0,1):
(a) for no nonempty open interval I in (0,1) is M_3 classical on I (Definition g8b:def:classical);
(b) M_3, continued holomorphically, has a length-four component S_4(M_3) that depends only on the function
    (Corollary g8b:cor:inv), equals the symbol S(Omega) of the table of sec_g8_M3alpha, and
    delta22(S_4(M_3)) = (1-1/c) p^r - (1+1/c) q^r != 0.
Supporting results: Lemma g8b:thm:indep (line 77): hyperlogarithms L^b_w (all words over a finite alphabet Sigma,
base point b in a simply connected domain U avoiding Sigma) are linearly independent over C(c); full proof
(minimal relation in lexicographic order (top length, number of top words), normalise, differentiate, residues).
Panzer (Comput. Phys. Commun. 188 (2015), Sec. 2.1, property 3 after Eq. (2.5)) cited as pointer only.

## Definition of the class C (line 9)
Classical factors on a simply connected domain U: (i) constants in C; (ii) branches l of log g (e^l = g),
g in C(c) non-constant, no zeros/poles in U; (iii) branches of Li_n(g), n>=1, defined as the last member of a chain
P_1..P_n of holomorphic functions with e^{-P_1} = 1-g, P_k' = P_{k-1} g'/g; g in C(c) non-constant, g and 1-g
without zeros/poles in U. C(U) = finite sums sum_i h_i prod_j A_ij, h_i in C(c). Any weight, any number of factors,
complex constants (zeta values etc.), arbitrary constants of integration in the chain (so every analytic continuation
of Li_n(g(c)) is included). "F classical on I": F = C on I (off poles of h_i) for some C in C(U), U a simply
connected domain containing I. The paper's Q(c)-combinations of products of Li_k(r), r a unit of Q(c), are included.

## Steps I consider delicate (new argument, needs referee)
1. Linear independence proof, sec_g8_M3func.tex:77-96. The minimality is over (N,k) lexicographic; the derived
   relation has g_u = 0 so it is trivial; then top coefficients constant; then residues at the distinct points of
   Sigma. Relations among words of different lengths are covered (the relation may contain any lengths).
2. Base point 0 is never used for L: H_w (base point 0, regular words) is written on any simply connected
   U in D\{0} as sum_{w=uv} s_u L^b_u(c) H_v(b) (Lemma g8b:lem:H, line 146); constants H_v(b) only multiply
   shorter words; s_u = (-1)^{#1s} from f_1 = 1/(1-t) = -1/(t-1). Base-point change is unitriangular with constant
   coefficients (Lemma g8b:lem:basic(b), line 57).
3. Ordering convention: rho(s_1...s_n) = e_{s_n} x ... x e_{s_1}, outermost letter last; S_{|w|}(H_w) = s_w rho(w)
   equals the section's g(a_n) x ... x g(a_1) with g(1) = -y, using e_1 = [c-1] = [1-c] mod constants (line 112ff, 146).
4. Constants: [g] = sum ord_s(g) e_s, kernel = constants; d log(const) = 0, so constant letters never occur
   (Lemma g8b:lem:rec, line 130). The section's "modulo the constant letter [2]" is exactly this: Lemma g8b:lem:two
   (line 241) shows the section's constant-free symbol of Li_{2,2}(1/2, 1+-c) IS the function invariant S_4,
   so the two-term identity of Lemma g8:prop:two holds for the (2,2) parts of the functions. The four-component
   remainder with [2] kept is not an invariant in this sense; nothing claimed.
5. C(c) vs Q(c) scalars: S_k is C(c)-linear, coefficients never differentiated in the invariant (they are
   differentiated only inside the independence proof). delta22 over K = C(c) via the field-agnostic proof of
   Lemma g8:lem:vanish; Q(c) -> C(c) and V_{Sigma_1} -> V_Sigma injective (line 173-175).
6. Length vs weight (line 170): constants have length 0; a classical C may have length > 4 (Li_7, products), so
   Lemma g8b:prop:classical (line 195) controls the length-four component at a FIXED base point, not the top one;
   this is what the comparison with M_3 needs. Length-k components of a branch of Li_n(g) (Lemma g8b:lem:classical,
   line 177): k=n: -[1-g] x [g]^{n-1}; k<n: P_{n-k}(b) [g]^{x k}; k>n: 0.
7. From an interval to a disc (proof of (a), line 212-232): choose c* in I off all exceptional points and a disc D'
   inside U, away from Sigma, inside D\{0}; identity theorem; enlarge alphabet (Corollary g8b:cor:unique(b)).
8. M_3 = R_3 - zeta(3)Phi + zeta(2)Psi: Phi = 1/2 - (1-c^2)/(4c)(H_1+H_{-1}) (length <=1, proved inline),
   Psi length <= 2, the 18 lower R_3 terms length <= 3. So S_4(M_3) = S(Omega) exactly; nothing discarded.
9. Bonus (line 234): Lemma g8:prop:one at function level, for G in L_Sigma (constant-coefficient hyperlogarithm
   combinations), plus any classical C.

## What is not claimed (line 40, gap box line 271)
Nothing at one angle; nothing for algebraic non-rational arguments (Li_n(sqrt c), Li_n(e^{i alpha})); not that
M_3 minus the two Li_{2,2} terms is classical (Goncharov's converse); no motivic identification of delta22;
nothing about the [2]-kept remainder; data of the 74/14 terms rest on Lemma g8:prop:R3 (M5).

## Machine check: scripts/g8b_check.py (1.6 s), output g8b_check.out
Part 1 exact: word lengths (R_3: 2/16/56 of lengths 2/3/4; Psi: const + 6 + 8 of lengths 1,2); rho(T_4) of the
56 terms equals the 24-row table of sec_g8_M3alpha (typed by hand) - True; delta22 formula - True.
Part 2 exact (sympy): S_4 of Li_n(g), n=3..6, five sample g, formal base-point constants: formula of
Lemma g8b:lem:classical True, delta22 = 0; six products (lengths 4..6) delta22 = 0.
Part 3 numeric (50 digits): expansions of Li_4, Li_5((1-c)/(1+c)) in L^{1/2}_w and Lemma g8b:lem:H for five words;
max difference 1.3e-51.

## Build
`tectonic -C --keep-logs wrap_g8b.tex` (inputs g6 flat/field/ladder/fourdilog, g8 M2alpha/M3alpha/M4alpha as on disk
at build time, then sec_g8_M3func). 0 undefined references, 0 overfull boxes, 0 "??" (pdftotext), no duplicate
hyperref destinations. 57 pages; the new section is Section 8, pp. 46-56.
Note for the merge: `supplement.tex` lacks a \bibitem{Panzer} (it is in refs.bib); wrap_g8b.tex adds one.
The main theorem is numbered "Theorem G8B" by a local \renewcommand{\thetheorem} inside a group.
