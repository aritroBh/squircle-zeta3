# COR8B — Explicit depth‑two part of M₃(α) at a general angle: settled

**Bottom line.** With the paper's own conventions (convB, δ₂,₂(a⊗b⊗c⊗d)=(a∧b)∧(c∧d), letters over (2,c,1−c,1+c)), COR8's two‑term answer is **exactly right, sign +, but only after dropping the constant letter [2]**. With [2] kept, DT_A (and DT_B — their δ₂,₂ are identical) misses by exactly four [2]-components, and I prove that **no** ℚ(c)-combination of depth‑two MPLs at {2,c,1±c}-unit arguments can repair it: the reachable δ₂,₂-space is only 3-dimensional and does not contain δ₂,₂(Ω). The obstruction is broken by allowing weight‑4 HPLs at S-unit arguments (then the target becomes reachable). COR8R's "−δ₂,₂(Ω)" was an artifact of computing with the other (convA) function.

---

## 1. Two‑variable symbols (Task 1) — computed and triply verified

Conventions (validated against series to 40 digits, `steps/step_01.py`):

Li_{2,2}(x,y)=Σ_{n>m≥1}xⁿyᵐ/(n²m²)=G_phys(0,1/x,0,1/(xy);1), Li_{3,1}(x,y)=G_phys(0,0,1/x,1/(xy);1).

- **Word identification (numeric)**: at (x,y)=(0.5,0.7), series vs iterated‑integral series agree to all 40 printed digits, e.g. Li₂,₂ = `0.06024453659644222914170370797481369677346`; shuffle Li₂,₂(x,y)+Li₂,₂(y,x)+Li₄(xy)=Li₂(x)Li₂(y) residual **−5.77·10⁻⁵³**; stuffle Li₃(x)Li₁(y) residual **−5.79·10⁻⁵³**.
- **dlog recursions (symbol backbone)** — each verified by finite differences at (0.5,0.7) and (0.31,1.4), residuals 10⁻¹⁸…10⁻²¹ (`steps/step_02.py`; note the correct dLi₁₁ = −Li₁(xy)[dlog y + dlog(1−x) − dlog(1−y)] − Li₁(x)dlog(1−y), which both COR8's referee and I initially got wrong):
  - dLi₂,₂ = Li₁,₂ dlog x + Li₂,₁ dlog y
  - dLi₃,₁ = Li₂,₁ dlog x − Li₃(xy) dlog(y/(y−1)) − Li₃(x) dlog(1−y)
  - dLi₁,₂ = −Li₂(xy) dlog(1−x) + Li₁,₁ dlog y
  - dLi₂,₁ = Li₁,₁ dlog x − Li₂(xy) dlog(y/(y−1)) − Li₂(x) dlog(1−y)
- **Symbols** (12 terms for Li₂,₂, 10 for Li₃,₁, 6 for Li₁,₃; saved in `twovar_symbols.json`): the two‑variable shuffle identities hold **exactly** (residual 0 after multilinear expansion).
- Engine validation end‑to‑end: S(F) computed from the 56 HPL terms of `hpl_fit_k3_(0,0,0).json` has **24 tensors** (= COR8R's count), is **byte‑identical to COR4's `F_symbol.json`** (0 differing tensors), and its δ₂,₂ = `(1−1/c)([c]∧[1−c])∧([1−c]∧[1+c]) + (−1−1/c)([c]∧[1+c])∧([1−c]∧[1+c])` — exactly the paper's Prop. prop:M3angle (2 components, printed coefficients {0:1,−1:−1} and {0:−1,−1:−1}).

General formulas (proved by the engine and used throughout): with X=[x],U=[1−x],Y=[y],V=[1−y],Z=[xy],W=[1−xy],
δ(Li₂,₂)=(W∧Z−U∧X)∧(Y∧V), δ(Li₃,₁)=−(W∧Z)∧(Y∧V)+(U∧X)∧(V∧X), δ(Li₁,₃)=(W∧Z)∧(X∧U). At x=½, δLi₃,₁ = −δLi₂,₂, so **δ(DT_A)=δ(DT_B) identically**.

## 2. Verdict on DT_A/DT_B (Task 2)

| Statement | Verdict | Evidence |
|---|---|---|
| δ₂,₂(Ω) = paper's formula | PROVED (re‑derived) | exact, residual 0 |
| One depth‑two term suffices | FALSE (as known) | coefficient ratio (c−1)/(c+1) nonconstant |
| DT_A = DT_B at δ₂,₂ level | PROVED | engine equality |
| [2] dropped: δ(DT_A)=δ₂,₂(Ω), sign **+** | **PROVED, exact** | residual **0** (COR8 right; COR8R's "−" came from convA) |
| [2] kept: δ(DT_A) = δ₂,₂(Ω) + R | PROVED, R has **exactly 4 [2]-components** | printed below |
| Any ℚ(c) combination of Li₂,₂/Li₃,₁/Li₁,₃ at S‑unit pairs matches δ₂,₂(Ω) with [2] kept | **IMPOSSIBLE (proved within the class)** | 3‑dim obstruction, functionals below |

The residual (`DT_A_residual_2kept.json`): R = (m∧p)∧(c∧p)·(c⁻¹−1) + (m∧p)∧(c∧q)·(c⁻¹+1) + (m∧q)∧(c∧p)·(−c⁻¹+1) + (m∧q)∧(c∧q)·(−c⁻¹−1) — precisely COR8's "4‑component [2]-residual".

**Obstruction (the new structural fact).** For all 344 S‑unit argument pairs (38 args, exhaustive within the Mason–Stothers degree bound; 720 terms) the δ₂,₂'s span only the 3‑dimensional space ⟨K∧([c]∧[1−c]), K∧([c]∧[1+c]), ([c]∧[1−c])∧([c]∧[1+c])⟩ with K = [2/(1+c)]∧[(1−c)/(1+c)]. Every term satisfies the 4 linear constraints (mp^cp)+(mq^cp)=0, (mp^cq)+(mq^cq)=0, (mp^cp)+(cp^pq)=0, (mp^cq)+(cq^pq)=0 (annihilator printed). The target violates the last two (evaluations f(T⁰)=1, f(T¹)=−1 and f(T⁰)=−1, f(T¹)=−1). Robustness: 1164 randomized large‑exponent S‑unit pair tests, **0 violations**. Hence the minimal exact set with constant letters dropped is **two terms** (DT_A or DT_B); with [2] kept, depth‑two MPLs at S‑unit arguments **never** suffice.

**Escape hatch found:** adding weight‑4 HPLs at the 9 S‑unit args with [1+g] in the alphabet, the target becomes reachable (combined span rank 7; "T0 in span: True, T1 in span: True" printed). So the [2]-residual is repairable by HPL terms, not by MPLs.

## 3. Function identity (Task 3)

- **[2]-dropped symbol identity PROVED and explicit** (confirms COR8's claim D, which COR8R left unconfirmed): S(Ω−DT_A) = Σ a_B(c) S(B) over Li₄'s and shuffle products; solved at c⁻¹, c⁰, c¹ with **reconstruction residual exactly 0** at all three exponents (68/66/72 terms; `classical_identity.txt`; my earlier projection bug — dropping vs projecting [2] — identified and fixed; after the fix COR8's own solution matches my target too).
- **The naive function identity Ω = DT_A + (Li₄'s+products) provably cannot hold**: δ₂,₂(Ω−DT_A−C)=R≠0 while every classical C has δ=0. Numerically at c=½ (`steps/step_11.py`): Ω = `−0.0416576932572766672882308026250919127289839427`, DT_A = `−0.370468167158478130188006668187879050052512399` (Ω−DT_A = 0.32881047390120144), while the principal‑branch classical combination gives 334.7−533.1i (mine) / 39.2+6.5i (COR8's — this is exactly COR8's unexplained "−38.9 gap", now explained: it is the [2]-residual R plus branch effects, not merely periods).
- A **genuine** function identity therefore needs the [2]-letter HPL terms (reachable per §2). The search for the explicit HPL combination X with δ(X)=−R (subset search over the 192 HPL δ's) **timed out** and is the only unfinished piece; the final PSLQ at two precisions depends on it and was not run. Recipe is in place: all evaluators (HPL series à la COR4, convB MPL series) are written and validated (`steps/step_11.py`).

## 4. Sentence the paper could carry

> "Modulo constant letters, the (2,2) obstruction of M₃(α)'s weight‑4 part is accounted for by exactly two depth‑two functions: δ₂,₂(Ω) = δ₂,₂( ((1−c)/c)Li₂,₂(½,1−c) − ((1+c)/c)Li₂,₂(½,1+c) ) (equivalently with +((1+c)/c)Li₃,₁(½,1+c)); the constant‑letter residual (m∧p)∧(c∧p)(c⁻¹−1)+(m∧p)∧(c∧q)(c⁻¹+1)+(m∧q)∧(c∧p)(1−c⁻¹)+(m∧q)∧(c∧q)(1+c⁻¹) is intrinsic — no depth‑two multiple polylogarithms at {2,c,1±c}–unit arguments reproduce it, since their δ₂,₂'s span only the 3‑dimensional space ⟨[2/(1+c)]∧[(1−c)/(1+c)]∧([c]∧[1±c]), ([c]∧[1−c])∧([c]∧[1+c])⟩ — but it is reproduced by weight‑4 HPLs at S‑unit arguments."

## 5. Files (all in `new_results/COR8B_depth_two_two_terms/`)

`cor8b_sym.py` (validated exact symbol engine), `twovar_symbols.json`, `DT_A_residual_2kept.json`, `pairs_delta_all.pkl`/`pairs_delta_all_v2.pkl` (MPL δ database), `hpl_delta_all.pkl` (HPL δ database), `syms_omega_dt.pkl`, `targets81_fixed.pkl`, `classical_sols_tracked.pkl`, `classical_identity.txt` (explicit [2]-dropped identity), `steps/step_01.py` (word ID, 40 digits), `step_02.py` (dlog checks), `step_06.py`, `step_07.py` (S‑unit δ databases), `step_08.py` (obstruction robustness), `step_11.py` (function‑level evaluation), `progress.log` (one line per step). Key evidence numbers are quoted above; nothing is claimed that was not printed.