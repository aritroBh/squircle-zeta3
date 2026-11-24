# RST1 — the rest exponent $\sqrt{2\pi}$ in the decay law of the squircle coefficients $M_k$

## Theorem (PROVED, computer-assisted)

For every integer **$N=k+1\ge 100$**,
$$
\boxed{\;|M_k-M_k^{\rm main}|\;\le\;1030\;N^{-1/2}\,e^{-\sqrt{2\pi N}}\;}
$$
(arb-certified raw constant $C_{\rm raw}=1025.889838913978\ldots<1030=C$; margin $4.110$, certified). Hence $|M_k-M_k^{\rm main}|\le 1030\,e^{-\sqrt{2\pi N}}$, i.e. the rest of Observation `obs:Mk` is $O(e^{-\sqrt{2\pi N}})$ with explicit constant $1030$ — the **true rest exponent** $\sqrt{2\pi}=2.5066282746\ldots$, previously only numerical (CX6 proved the weaker $e^{-2\sqrt N}$; $2<\sqrt{2\pi}$). This closes the last numerical piece of the decay law.

**Corollary (full $M_k$, improving DCY2).** With DCY2's main term ($C_2=10^6$, $\varpi_1$):
$$\Bigl|M_k+\tfrac{\sqrt\pi}{N}e^{-\theta}\operatorname{Re}\bigl[e^{i\theta}\bigl(1+\tfrac{\varpi_1}{\sqrt N}\bigr)\bigr]\Bigr|\le C_3 N^{-2}e^{-\theta},\qquad C_3=1\,000\,668\quad(N\ge100),$$
since the new rest contributes $\le 667.357\,N^{-2}e^{-\theta}$ ($\max_{N\ge100}1030N^{3/2}e^{-(\sqrt{2\pi}-\sqrt\pi)\sqrt N}=667.357$ at $N=100$, decreasing for $N>16.7$; arb + grid). DCY2 had $C_3=1.85\times10^{10}$ from CX6's rest.

## Item 1: which singularity produces the rest

**The rest is produced by the pole of the weight kernel at $s=-2$ — the $m=2$ image of the same $\pm i\pi/2$ Borel branch points — NOT by the next branch points $\pm3i\pi/2$.** In CX6's representation, the main term pairs the image $E=e^{i\pi s/2}$ ($s=1$ pole of $K$; $z=i\pi/2$) with the $Q_k$ pole at $s=-1$ (decay $e^{-\sqrt{\pi N}}$). The rest pairs $E$ with the **next kernel pole $s=-2$** through the explicit $q=2$ piece $A_k^{[2]}(s)=-\frac{2s^{k-1}(s+1)}{(s+2)^{k+1}}$ of CX6's $A_k$: saddle $r=2\sqrt{N/\pi}$, decay $e^{-\sqrt{2\pi N}}$. The next branch points ($s=\pm3$, $z=\pm3i\pi/2$, the $K-E\sim-E^3$ term) and the $s=-3$ pole ($m=3$) each give only $e^{-\sqrt{3\pi N}}$ ($\sqrt{3\pi}=3.0699>\sqrt{2\pi}$).

**Leading asymptotics (saddle of the $q=2$ model, as DCY2 for the main pair):**
$$
M_k-M_k^{\rm main}\sim R_k:=-\frac4\pi\operatorname{Im}\int_{\rm ray}H(s)e^{i\pi s/2}\frac{2s^{k-1}(s+1)}{(s+2)^{k+1}}\,ds
\;\sim\;\frac{\sqrt\pi}{N}e^{-\sqrt{2\pi N}}\cos\sqrt{2\pi N}.
$$

## Lemma ledger

| # | Lemma | Status | Key number |
|---|---|---|---|
| L0 | CX6 inputs: decomposition (3); $|A_k|$ bound; $|H|,|K|,|K-E|$ bounds; CX6 (17),(18) | **PROVED** (CX6, refereed; reused) | — |
| R1 | Modulus-exponent bound $\tilde E(r)=-ar-\frac{\sqrt2 jq}{2r+c_0}$ | **PROVED (hand)** | Padé $\log(1+u)\ge\frac{2u}{2+u}$, AM–GM |
| R2 | Gaussian $F_{H,\lambda}$ saddle lemma ($r_\lambda>c_0$) | **PROVED (hand)** + **CERTIFIED** | exact/bound $\le 0.00165$ |
| R3 | Universal box lemma | **PROVED (hand)** + **CERTIFIED** | exact/bound $\le 0.00499$ |
| R4 | $T_1$ assembly ($q2$/mid$[2,12]$/mid2$[13,k/2]$/tail$>k/2$) | **PROVED (hand)** + **CERTIFIED** | — |
| R5 | Convexity $Q_i\le Q_2$ on $[2,k/2]$, $Q_i\le Q_{13}$ on $[13,k/2]$ | **PROVED (hand)** + **CERTIFIED** | $h(99)=12.78>0$; 0 violations |
| R6/R7 | $T_2$ (decay $e^{-\sqrt{3\pi k}}$); $T_3,J_k$ elementary | **PROVED** | — |
| R8 | Every normalized piece non-increasing for $k\ge99$ | **PROVED (hand)** + **CERTIFIED** | 0 grid violations |
| R9 | Constant $C_{\rm raw}<1030$ | **CERTIFIED (arb)** | $1025.8898389140$ |
| R10 | Sanity $k\in[99,400]$ ($N\in[100,401]$) | **CERTIFIED (arb quadrature)** | max ratio $0.215$ |
| R11 | Item 1: rest $=R_k(1+o(1))$; constants; subdominance | **CERTIFIED (numerics)** | below |

## Evidence table (printed numbers)

**Item 1 (step_03, step_19, step_06).** rest$/R_k=0.99971,1.00262,1.00004,0.99981,0.99999,0.99997,0.99999$ at $k=100,150,200,250,300,350,400$ (exact $M_k$ from `mk_deg9.tsv`). Saddle amplitude $|g_k|Ne^{\theta_2}=1.98950,1.86506,1.81302,1.79043$ at $k=400,2000,10^4,5\times10^4$ ($\to\sqrt\pi=1.7724539$); phase $\arg(g_k)-\theta_2=1.51096,1.54405,1.55884,1.56545$ ($\to\pi/2=1.5707963$). Subdominance: $|R^{(q=3)}|/|R|=9.8\mathrm{e}{-4},4.8\mathrm{e}{-5},4.0\mathrm{e}{-6}$ and $|R^{(K-E)}|/|R|=4.5\mathrm{e}{-3},2.5\mathrm{e}{-4},3.4\mathrm{e}{-5}$ at $k=100,200,400$.

**Constant pieces at $k=99$ (step_15 / verify_rst1 Part A; arb, 80 digits).**

| piece | value (arb) |
|---|---|
| $q=2$ ($i=0,1$) | $1012.6418023331$ |
| mid ($i\in[2,12]$) | $2.6180577169686$ |
| mid2 ($i\in[13,k/2]$) | $5.6494647\times10^{-17}$ |
| tail ($i>k/2$) | $1.8397458\times10^{-5}$ |
| $T_2$ | $10.629960466490$ |
| $T_3+J_k$ | $1.2138883\times10^{-18}$ |
| **$C_{\rm raw}$** | $1025.8898389140$ |

$a=1.1107207345395916$, $\sqrt{2\pi}=2.5066282746310005$, margin $4.1101611$.

**Lemma validation (step_18b / Part B; exact $F_{H,\lambda}$ by arb quadrature).** worst exact/bound ratio $=0.00165$ (Gaussian R2), $0.00499$ (box R3) over $(j,q,\lambda)$ including $(99,2,a),(98,2,a),(400,2,a),(87,13,a),(97,3,a),(91,9,a),(99,1,3a),(400,1,3a)$. Convexity: $Q_i\le Q_2$ on $[2,k/2]$ (0 violations), $Q_{13}\ge Q_{k/2}$ on $[99,400]$ (0), $h(99)=12.782895>0$, $h(k)>0$ for $k\ge99$. Monotonicity (Part D): 0 violations on $[99,1000)$.

**Certified sanity table (step_17 / Part E; arb, 302 quadratures).** All $k\in[99,400]$: $|M_k-M_k^{\rm main}|<1030N^{-1/2}e^{-\sqrt{2\pi N}}$; **max $|{\rm rest}|/{\rm target}=0.21501950$ at $k=106$**.

| $k$ | $|M_k-M_k^{\rm main}|$ (arb) | /target | /lead ($\frac{\sqrt\pi}{N}e^{-\theta_2}$) |
|---|---|---|---|
| 99 | $1.64857408084\mathrm{e}{-13}$ | 0.127 | 0.716 |
| 150 | $4.403845828\mathrm{e}{-17}$ | 0.013 | 0.089 |
| 200 | $3.811407939\mathrm{e}{-18}$ | 0.147 | 1.174 |
| 250 | $1.884388565\mathrm{e}{-20}$ | 0.053 | 0.471 |
| 300 | $1.927871694\mathrm{e}{-22}$ | 0.026 | 0.252 |
| 350 | $1.220912698\mathrm{e}{-23}$ | 0.057 | 0.601 |
| 400 | $4.799125366\mathrm{e}{-25}$ | 0.061 | 0.684 |

(/lead oscillates in $[0.09,1.18]$ as $\cos(\theta_2+\varphi_N)$ — the predicted ratio.)

## Proof route (as tasked)

Move the contour past the next singularities: CX6's exact decomposition isolates the rest; the **next pair** ($\pm3i\pi/2$, the $K-E$ term $T_2$) is bounded explicitly with sharp exponent $e^{-\sqrt{3\pi N}}$ via the Gaussian $F_{H,3a}$ lemma (as DCY2 bounded the main pair); the **rest-producer** ($s=-2$ pole, the $q=2$ piece of $A_k$) is bounded explicitly with the sharp exponent $e^{-\sqrt{2\pi N}}$ via the Gaussian $F_{H,a}$ lemma (power $N^{-1/2}$ from the modulus bound, no phase cancellation); everything else is bounded by $e^{-c\sqrt N}$ with $c>\sqrt{2\pi}$ (mid2: $\sqrt{14\pi}=6.63$; tail: $e^{-k/2}$; mid $[2,12]$: $\sqrt{3\pi}$). Each normalized piece is non-increasing for $k\ge99$, so the constant is the sum of the arb-evaluated pieces at $k=99$.

## Files (work dir `new_results/RST1_rest_exponent/`)

- `proof.md` — full proof with the lemma ledger (L0, R1–R11) and all derivations. `result.md` — this report.
- `verify_rst1.py` — master verifier: recomputes every certified constant (arb), validates R2/R3 vs exact arb quadrature, checks convexity and monotonicity, rebuilds $C$, re-runs the sanity table $k\in[99,400]$, re-checks item 1. **Full run: PASS, 282 s.**
- `step_03.py`–`step_06.py` — discovery numerics (rest law, saddle constants, phase). `step_07.py`–`step_14.py` — $F$-lemma development + grid validation of bounds/monotonicity. `step_15.py` — arb constants ($C_{\rm raw}$). `step_16.py`, `step_17.py` — certified sanity table (checkpoint `sanity_progress.json`). `step_18.py`, `step_18b.py` — certified validation of R2/R3 (exact arb quadrature) + convexity. `step_19.py` — item 1.
- `rst1_common.py`, `rst1_common2.py` — shared numerics; `progress.log` — step log.
- Inputs (read-only): `D13f_sweep/mk_deg9.tsv` (certified $M_k$ balls), CX6's lemmas (`CX6_decay_rest/proof.md`), DCY2 quadrature settings (`DCY2_decay_first_correction/step_06_finite.py`).

Reproduce: `cd new_results/RST1_rest_exponent && '<project>/.venv/bin/python' verify_rst1.py` (PASS, ~5 min).