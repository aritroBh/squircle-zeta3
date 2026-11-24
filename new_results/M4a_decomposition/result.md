# M4a — weight‑5 (W5) and weight‑4 (W4″) parts of M4: ALL 5 TASKS COMPLETE

**Verdict: PASS on every task.** (1) M4 = W5 + W4″ + W3⁽⁴⁾ + T verified to 5.5×10⁻¹⁰⁸, with M4 independently confirmed by Neville extrapolation (4.4×10⁻⁵⁶). (2) Exact G‑expressions for W5 (782 terms) and W4″ (432 terms) written as TSVs, each matching quadrature to ~200 digits. (3) Level‑8 word vectors built exactly (W4″: 3524 rational + 2638 √2 words; W5: 18900 + 25762), pickled, and matching quadrature to ~10⁻⁶². (4) W4″ to 800 digits (two‑run agreement 2.7×10⁻⁸¹⁰), W5 to 600 digits (8.3×10⁻⁶¹²), M4 to 800 digits (3.6×10⁻⁸¹³), plus W4″ to 1060 digits. (5) A candidate closed form for W4″ in the same 60‑number basis as W4, stable at 700 and 1000 digits (identical up to global sign; residuals 6.9×10⁻⁷¹⁶ and 1.5×10⁻¹⁰¹⁶ vs random sizes 10¹¹·⁵ and 10¹⁶·⁴), agreeing with the 1060‑digit W4″ to 1.07×10⁻¹⁰⁶⁰.

## Important operational note (false alarm resolved)
`import dilog_exact` (triggered by unpickling the paper's `*_tails_exact.pkl`) executes `mp.mp.dps = 60` at import time, so any mpmath evaluation after `pickle.load` silently runs at 60 digits. Every apparent ~10⁻⁶⁰ "discrepancy" with the paper's exact values was this artifact. With dps set after imports: the paper's reduced W3′ (k=3) and W3⁽⁴⁾ (k=4) agree with quadrature to 1.2×10⁻¹¹⁰; the rem:k4 identity W3⁽⁴⁾ = W3′(k=3) − ζ(3)/8 − π²λ/24 − λ³/12 holds to 2.8×10⁻¹²⁶⁴ (flint prec 4200); and M3 − W3′ − T equals the proved thm:W4 closed form to 2.2×10⁻¹¹⁹⁸ (all‑flint recomputation). **The paper's values are correct.**

## 1. Exact Ω's for k = 4 (re‑ran `ibp_reduce.py 4`; output byte‑identical to paper's `ibp_reduce_k4.out`; parsed to pole/residue data `omega4.json` with exact sympy reconstruction asserted)

With a = 1+√2, b = √2−1, M4 = Σ ∫_b^a l1ⁱ l2ʲ Ω_ij(v) dv:

- **W5** (i+j=4): Ω₂₂ = (7√2/128)/v; Ω₃₁ = (−5√2/96)/v; Ω₄₀ = (1/24 − 3√2/128)/v + (−1/12)/(v+a).
- **W4″** (i+j=3): Ω₂₁ = (7/16+√2/8)/v + (5/16)/(v−a) + (−7/8)/(v+a) + (−5/16)/(v+b); Ω₃₀ = (1/24−√2/8)/v + (−1/8)/(v−a) + (5/24)/(v+a) + (−1/6)/(v+b).
- W3⁽⁴⁾ (i+j=2): Ω₁₁ = (3/4+√2/4)/v + (−3/2)/(v+a); Ω₂₀ = (1/4−√2/4)/v + (1/4)/(v−a) + (1/4)/(v+a) + (−1)/(v+b).
- T (i+j=1): Ω₁₀ = (1/2)/v + (1/2)/(v−a) + (−1)/(v+a) + (−1/2)/(v+b).

## 2. G‑expressions (files `gpl_w5.tsv`, `gpl_w4pp.tsv`, same format as `A3_W4_proof/gpl_expression.tsv`)
Built with `w4_symbol.log_product` from (pole z, l1ⁱl2ʲ factors, Q(√2) coefficient) triples parsed from `omega4.json`; every entry satisfies len(word)+k = weight; alphabet {0,±a,±b}. Checks: own copy of the `real_g` acb Taylor evaluator **plus shuffle endpoint‑regularization** (unit tangent G_reg((e,);e)=0 — needed because W4″ has poles at z=a; 108 such words), prec 700 bits, cap 420:

| piece | terms | \|G − quad\| (90 quad digits) | flint radius |
|---|---|---|---|
| W4″ | 432 | 0.0 | 2.65×10⁻²⁰⁸ |
| W5 | 782 | 0.0 | 7.99×10⁻²⁰⁹ |

## 3. Level‑8 word vectors (files `wv_w4pp.pkl`, `wv_w5.pkl`)
Converted exactly as `verify_D.rebuild_target`: `forms_path` from POINTS[1] to POINTS[A] minus POINTS[B], shuffled with `atom(('log',2))^k`, split into rational and √2 parts. All words have correct weight. Numeric check (own copy of `check_conversion`'s series/half/g_at_one/h evaluator, prec 400, cap 330):

| piece | support (rat, √2) | \|wv − quad\| | flint radius |
|---|---|---|---|
| W4″ | (3524, 2638) | 9.52×10⁻⁶³ | 7.8×10⁻¹¹⁶ |
| W5 | (18900, 25762) | 8.14×10⁻⁶² | 1.6×10⁻¹¹⁵ |

## 4. High‑precision values (mpmath tanh‑sinh on the v‑form = mk_numerics method A; two independent working precisions; timing extrapolated from a 300‑digit probe)

| value | runs | \|diff\| | file | digits |
|---|---|---|---|---|
| W4″ | @890/@850 | 2.67×10⁻⁸¹⁰ | `W4pp_800.txt` | 800 |
| W5 | @690/@650 | 8.30×10⁻⁶¹² | `W5_600.txt` | 600 |
| M4 | @890/@850 | 3.63×10⁻⁸¹³ | `M4_800.txt` | 800 |
| W4″ | @1100 (err est 10⁻¹¹⁶⁴) | — | `W4pp_1060.txt` | 1060 |

W4″ = 1.424453676008872921514422491257270535637099847675524471725467791684571969152806…
W5 = −0.18235662466773630012609370629766450393982143738931837189657087268442643560389849…
M4 = 0.006027993112609122504516108250877605571476734148292802895459093539190168935771…
Independent‑method evidence: G‑expression values (Task 2), word‑vector values (Task 3), Neville M4 (4.4×10⁻⁵⁶), final sum W5+W4″+W3⁽⁴⁾+T−M4 = −2.8×10⁻⁶⁰¹ (limited by W5's 600‑digit truncation).

## 5. Candidate closed form for W4″ (λ = log 2, L = log(1+√2), b = √2−1)
Same 60‑number basis as the W4 search (22 weight‑4 products ×2 plus the same 8 Li4 arguments ×2). LLL relation at 700 and 1000 digits **identical up to sign** (W4″ coefficient 11520, max|coef| 1.198×10⁶, residuals −6.87×10⁻⁷¹⁶ / −1.45×10⁻¹⁰¹⁶; random sizes ~10¹¹·⁵ / 10¹⁶·⁴). A second, auto‑greedy‑basis 700‑digit run gave a consistent different‑basis relation (W4″ coef −6912, residual 3.43×10⁻⁷¹⁷). Verified: closed form − W4″(1060‑digit quadrature) = **1.07×10⁻¹⁰⁶⁰**.

W4″ = (1−104√2)Li₄(b) − (1+104√2)Li₄(−b) + 2√2 Li₄(1/√2) + (41/4)√2 Li₄(−b²) − (37/4+√2/8)Li₄(1/2) + (10+16√2)Li₄(1−1/√2) − (10−16√2)Li₄(2−√2)
+ L[(10+2√2)Li₃(b) − (10−2√2)Li₃(−b) − 18Li₃(1/√2) + (63/32+3√2/8)ζ(3)]
+ λ[6√2(Li₃(b)−Li₃(−b)) − 11√2 Li₃(1/√2) + (1/16+77√2/64)ζ(3)]
+ 5Li₂(b)² + Li₂(b)[(5+√2)L² − 4λL − (7/4)√2 λ² − (19/24+√2/6)π²]
+ (1333/11520 − 391√2/2880)π⁴
+ π²[(49/48+11√2/48)L² + (7√2/6−23/16)λL + (19/48−157√2/192)λ²]
+ (−19/24−√2/12)L⁴ + (−7/6−7√2/3)L³λ + (3/2+√2/8)L²λ² + (7/12+13√2/12)Lλ³ + (−37/96+13√2/48)λ⁴.

As expected, its π⁴ and Li/log coefficients differ from W4's while the basis is the same (no Li₄(1/√2−1) term, as for W4). This is a **candidate** (lattice‑found, numerically stable at two precisions, 1.07×10⁻¹⁰⁶⁰ agreement) — ready for the A3h5 exact certificate, whose word‑vector inputs are the pickles above.

## Step files that produced these results (all in `new_results/M4a_decomposition/`)
`step_quad_pieces.py` (Task 1 quadrature), `step_neville_m4.py` (method C), `gexpr.py` (Task 2 build), `g_eval.py` (Task 2 check), `step_wordvectors.py` (Task 3 build), `step_wv_check.py`/`step_wv_check5.py` (Task 3 checks), `step_hp_probe.py`/`step_hp.py` (Task 4), `step_pslq.py`/`step_pslq2.py` (Task 5 LLL), `step_assemble_w4pp.py` (closed‑form verification). Outputs: `ibp_reduce_k4.out`, `omega4.json`, `quad_90.json`, `m4_neville.json`, `gpl_w5.tsv`, `gpl_w4pp.tsv`, `*_check.json`, `wv_*.pkl`, `hp_results.json`, `W4pp_800.txt`, `W5_600.txt`, `M4_800.txt`, `W4pp_1060.txt`, `pslq*_w4pp_*.json`, `w4pp_closed_form.json`, `result.md`, `progress.log`. Rerun commands are listed in `result.md`.