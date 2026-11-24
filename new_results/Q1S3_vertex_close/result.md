# Q1S3: closed forms of V1 and V3' — status report

## Established (all printed by code)

### 1. Exact stage-2 reduction (PROVED algebra + NUMERICAL certification)
y=1/X maps X∈[1,∞) to y∈(0,1]. V1 = -G1 - b G2 + a G3, V3' = -H1 + (b√2/2)H2 + (a√2/2)H3 + √2 b L H4 + √2 a L H5
with y-form pieces (kernels 1/y, 1/(1±·y), 1/(1±·y)^2). Square kernels are exact derivatives:
G2,G3,H2,H3 are weight ≤2, H4,H5 weight 1; the weight-3 symbol of V1 (resp. V3') comes from G1 (resp. H1) alone.
Exact engine (q3_engine.py): everything reduces to Chen words G(w;1), letters {0,a,-b,±ia,±ib} (a=1+√2, b=√2-1),
plus boundary log-products. Every rule control-tested; every piece (G1,G2,G3,H1..H5) matches quadrature exactly.
CERTIFICATION vs 1250-digit references (q3_eval.py GL-panel evaluator, 128 words):
- V1: residual 2.97e-121 (120.5 digits)
- V3': residual exactly 0 (150 digits)
V1: 96 w3-words + 50 w2-words + 23 products. V3': 72+30+10.

### 2. Symbols (q3_symbol.py, Goncharov recursion; all letter differences S-units: 40 pairs, 0 failures)
- S(V3') = rational, conjugation-invariant, 16 terms: (a,a,a):4 (a,p2,a):16 (a,q1,a):6 (a,q2,a):6 (p2,a,a):16
  (p2,p2,a):32 (p2,q1,a):12 (p2,q2,a):12 (q1,a,a):16 (q1,a,p2):4 (q1,p2,a):28 (q1,q1,a):3 (q2,...) symmetric.
  EXACTLY equals the previous session's X-word symbol symX_V3p.pkl (independent presentation) — cross-validation.
- S(V1) = A + sqrt2 B; A (6 terms), B (18 terms), both conj-invariant; no imaginary parts.

### 3. Span test + exact solves (PROVED at symbol level)
rank[493 Li3]=40, rank[Decomp]=44, rank[all]=64; rank[272 ReLi3]=24, rank[ReLi3|Decomp]=55.
All three symbols in the span. Unique solution in a 55-column independent basis (q3_solve_min.pkl):
- S(V3'): 16 ReLi3(z) terms + 10 decomp terms (coeffs to 36); identity verified EXACTLY (max deviation 0).
- S(V1)_A: 8 ReLi3 + 3 decomp (coeffs to 74). S(V1)_B: 16 ReLi3 + 10 decomp (coeffs to 416).
Decomp generators realized as log u * X(v,w), X(v,w) = (1/2) log v log w + Li2-combination with S(X)=v x w
(exact check per pair; q3_Xli2.pkl). ReLi3 generator = S(z)+conj S(z) -> value 2 Re Li3(z).

### 4. Value-level identity for V3' (NUMERICAL, 400 digits)
V3' = sum c_i 2 ReLi3(z_i) + sum d_k log u_k X(v_k,w_k) + W2 + P + R,
W2 = 1.182296014243883... (weight-2 words), P = 0.313685414311214... (boundary products), both real.
R = 114.397510812934955838... - 80.549742763361728117... i   (symbol-0, PROVED; value at 390 digits).
Pipeline end-to-end tested on an artificial combination: remainder 0 exactly (80 digits).
OPEN: R's closed form not found — PSLQ at 380 digits against ~120-element bases
(zeta(3), L(chi8,3), pi^3, sqrt2 pi^3, Cl3(pi/4,3pi/4), pi*G, pi*Cl2(pi/4), zeta(2)*(lam,L,l3,alpha),
pi*(lam^2,L^2,lam L,l3^2,lam l3,L l3,alpha^2,lam al,L al,l3 al), pi*(Re/Im)Li2 at alphabet + at 46 S-unit
cross-ratios, pi*Ti2(sqrt2), pi*Cl2(2 alpha), (Re/Im)Li3 at the 7 torsion-pair points) returned nothing.
The large coefficients of the current (pivot/minimal) ReLi3 combination (up to 36) make R huge
(cancellation -122.35 -> 2.08); a kernel-reduced short combination should shrink R and make the PSLQ tractable.

### 5. NOT done
Invisible remainder R for V3' (open); V1 realization (symbol solved, value remainder not computed);
short combination via kernel LLL; Q1 and s3 assembly; verify_q1s3.py.

## Files
q1s_exact.py, q1s_symbol.py (exact F + symbol engine, copied from Q1S_vertex_symbol),
q3_engine.py (stage-2 reduction), q3_eval.py (GL-panel Chen-word evaluator), q3_stage2.py (piece assembly),
q3_symbol.py (symbols of Acc), q3_accs.pkl, q3_symbols.pkl, q3_solve.pkl, q3_solve_min.pkl, q3_sel.pkl,
q3_Xli2.pkl, K_V3p_400.pkl, K2_V3p_400.pkl, R_V3p_400.pkl, sunit_li3.pkl.
