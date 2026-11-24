# DCYA1 — RESULT

**Goal:** prove $|M_k(\alpha)+(\sqrt{2\alpha\sin\alpha}/N)e^{-\sqrt{2\alpha N}}\cos\sqrt{2\alpha N}|\le C(\alpha)N^{-3/2}e^{-\sqrt{2\alpha N}}$.

**Status: partial — the full proof architecture is established and validated numerically; two explicit constants are not yet arb-certified.**

## PROVED
- **Exact representation** $M(u;\alpha)=T(u)+R(u)$ (generalises CX6's kernel representation to every angle), derived by hand, verified to ~50 digits at 8 $(u,\alpha)$ points. Branch points of $\widetilde B_\alpha$ located exactly at $z=\pm i\alpha+2\pi i n$ (source of the rate $\sqrt{2\alpha}$).
- **Weight expansion** $M_k=\sum_{m\ge1}(-1)^{m-1}T_k^{(m)}+R_k$.
- **Main term** reduced to explicit 1-D integrals $I_+^{(j)}$ (Tricomi-$U$); leading $+g_1$ captures the first correction $\kappa_1(\alpha)$ exactly (verified numerically: remainder $\sim N^{-2}e^{-\theta}$).

## NUMERICAL (validated, not arb-certified)
- Saddle model $MV$ (linear phase $\omega$) approximates the main term to $O(N^{-3/2}e^{-\theta})$ ($|I_0-MV|/{\rm sc}\approx1.7$); $MV\to$ lead; saddle series $F_1,F_2$ capture the remainder.
- Decay-law residual ratio is $O(1)$ at $\alpha=\pi/3,2\pi/3$ ($k\le100$); $\alpha=\pi/2$ matches D13/K2.

## OPEN (gap, stated exactly)
- **L6:** arb-certified main-term constant $C_1(\alpha)$ (D13s-style subinterval cell sweep of the derived error function $E(N)$) — not completed.
- **L7:** arb-certified rest constant $C_2(\alpha)$ with $\delta=\sqrt{2(2\pi-\alpha)}$ (CX6-style domination at general $\alpha$) — machinery derived, constant not assembled.
- Uniform-in-$\alpha$ version and residual checks at $k=200,300$ — not reached.

## Files
proof.md (lemma ledger + derivation), verify_dcya1.py, mk_eval.py (reliable $M_k$ evaluator), mk_*.json ($M_k$ tables), progress.log, step_*.py.
