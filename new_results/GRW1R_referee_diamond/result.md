Adversarial referee report on GRW1 (diamond end). Full write-up in `result.md`; master script `verify_grw1r.py` (52 s) reproduces every number below; one line per step in `progress.log`; step files `step_01.py … step_36.py` (the main computational steps; a few exploratory/debugging files in between carry the numbering gaps — the canonical pipeline is cited per claim).

## 1. Verdict table

| Claim | Statement | Verdict |
|---|---|---|
| **A** | Laplace form $I(\eta)=\int_0^\infty e^{-u}H(\eta u)K(\eta,u)\,du$, $\eta>-1$ | **PROVED** |
| **B** | $f_n=\sum_{j=0}^n H_{n-j}M_{n,j}$ exactly; $M_{n,0}=n!\,\eta_D(n)$; $\lvert M_{n,j}\rvert\le (n-j)!/(j!(j+1)^{n-j+1})$ | **PROVED** |
| **C** | Darboux transfer for $H$: $E^*=5/24$ on the slits; $\lvert R_n\rvert\le 4.45\rho^{-n}n^{-7/2}+135(2.9\rho)^{-n}$ | **PROVED AFTER REPAIRS** — $E^*=5/24$ now certified (arb cover + exact series, replacing the 4000-point grid); phase typo in the printed Hankel display fixed; remainder derivation defective (false $1.63$ intermediate at $n=10$; unjustified absorption), repaired constant $7.06$; the stated $4.45$-bound is *true* (verified exactly $n\le240$, worst 0.986 at $n=18$) but not proved by the printed argument |
| **D** | Uniformity: $\lvert f_n-H_nn!\eta_D(n)\rvert\le C_He^{2\rho}\rho^{-n}n!2^{-n-1}/n$ | **PROVED** (given C as repaired) |
| **E** | Growth law $f_n=\sqrt2\cos\frac{\pi(3-2n)}4\,n!(2/\pi)^nn^{-3/2}\bigl(1+\frac{3-3\pi(-1)^n}{8n}+O(n^{-2})\bigr)$; rhombus law with $\theta_0$; radius zero | **PROVED AFTER REPAIRS** (inherits C's repairs; the $O(n^{-2})$ term *is* proved — relative constant $\le17.2$ for $n\ge10$) |
| **F** | Borel summability along $\arg\zeta=0,\pi$; $S(\eta)=I(\eta)$ for every $\eta>-1$ | **PROVED AFTER REPAIRS** — Lemma 6's "$J_0$-series strip extension" is **false** (the series diverges off $\mathbb R$); repaired via Lemma 6b (whose contour arithmetic must be written out); two cosmetic fixes |
| **G** | Singular directions; decomposition of $B$ in $\lvert\operatorname{Im}\zeta\rvert<\pi$ | **PROVED** with GRW1's own caveat (possible $R_2$ line-boundaries on $\lvert\operatorname{Im}\zeta\rvert=(j+1)\pi/2$, $j\ge1$, not pinned down; "directions exactly $\pm\pi/2$" holds modulo that; $\arg\zeta=0,\pi$ are *proved* nonsingular — what the paper uses) |
| **H** | Stokes jump $2\sqrt\pi\,e^{-3i\pi/4}\varepsilon^{1/2}e^{-\pi/(2\varepsilon)}(1-3i\varepsilon/4+O(\varepsilon^2))$ | **PROVED** (given the proved singularity structure); asymptotics verified, next coefficient $-25/32$ extracted |

## 2. The repairs (written out)

**R1 (C).** *$E^*=5/24$ is now a theorem* (was: 4000-point grid). $\psi(s)=\lvert E(is)\rvert^2=\Psi(s^2)$ with $\Psi(u)=\sum c_mu^m$ computed exactly ($c_0=(5/24)^2$, $c_1=-25/13824<0$, $c_2=11147/232243200$). Piece A $[0,0.1]$: $\Psi(u)\le(5/24)^2+(c_1+(c_2+T)u)u$, $T=M_E(3)^2/(9^3(1-u/9))$, $M_E(3)\le(\sqrt{(1+e^6)/6}+1+1.5)/9\le1.1900048$; arb: $c_1+(c_2+T)(0.01)\le-1.7885\mathrm e{-3}<0$. Piece B $[0.1,1.5]$: 45-term exact series + Cauchy tail on 700 arb intervals: $\max\psi\le0.0433847019$. Piece C $[1.5,1.9\rho]$: closed form on 744 arb intervals: $\max\psi\le0.0399729017$. Both $<(5/24)^2=0.0434027778$. **$\sup\lvert E\rvert=5/24$, attained only at $s=0$.**

**R2 (C).** *Phase typo in the printed Hankel display.* With GRW1's literal $2\operatorname{Re}[e^{i\pi/4}i^{-n}]$ phases the main term has the wrong sign pattern at even $n$. Correct phases (verified $\lvert H_n-A_n-B_n\rvert\le$ $E$-term bound at every $n$ tested, odd and even):
$$A_n=-\tfrac{2}{\pi}\cdot 2\cos(n\pi/2+\pi/4)\,\rho^{1/2-n}B(\tfrac32,n-\tfrac12),\qquad B_n=\tfrac{2}{\pi}\cos(n\pi/2-\pi/4)\,\rho^{3/2-n}B(\tfrac52,n-\tfrac32).$$
GRW1's final two-term law is exactly what these give — a writeup typo, not a math error.

**R3 (C).** *Remainder derivation defective; repaired.* The printed "$\Gamma(n-\frac52)/\Gamma(n+1)\le1.63\,n^{-7/2}$" is false at $n=10$ (true ratio $1.6306839$); the claim that the Beta $O(n^{-2})$ truncations are "absorbed in the same $4.45$" is wrong as accounting ($\lvert A_n+B_n-\mathrm{law2}\rvert\le2.7882\,u$ at $n=11$ is separate and additive). Honest assembly, arb-certified at every integer $10\le n\le500$ (sup $7.0523$ at $n=11$) + rigorous Stirling tail ($\le5.19$ for $n\ge500$):
$$\lvert R_n\rvert\le 7.06\,\rho^{-n}n^{-7/2}+135(2.9\rho)^{-n}\quad\text{for all }n\ge10.$$
The *stated* $4.45$-bound is true (exact rational $H_n$, worst ratio $0.9864$ at $n=18$ on $10\le n\le240$; asymptotic constant $\to3.98$) but not proved by the printed argument. Downstream constants $C_H=5.874$, $0.968\le\lvert H_n\rvert/\mathrm{lawH}\le1.032$ ($n\ge190$) are unaffected.

**R4 (F).** *Lemma 6's "strip extension" is false*: for $w$ off $\mathbb R$, $\lvert J_0(2\sqrt{(1-2k)w})\rvert\sim e^{c\sqrt k}$ ($c>0$) against $\lvert c_k\rvert\sim k^{-3/2}$ — the $J_0$ series diverges for every non-real $w$; it represents $H^{[0]}$ only on $\mathbb R$. Repair: strip/sector bounds come from **Lemma 6b**; in Lemma 8 Step 3 read "strip form of Lemma 6b". Lemma 6b itself needs its contour arithmetic written out: the representation $H^{[0]}(w)=\frac1{2\pi i}\oint e^sH(w/s)\,ds/s$ starts from a **large** circle $\lvert s\rvert>2\lvert w\rvert/\pi$ ($s=0$ is an *essential* singularity), deformed to the rectangle hugging the segment $\{tu\}$, $T=2\lvert w\rvert/\pi+1$, sides at $d=\sqrt{\lvert w\rvert}$. The argument is sound; the conclusion $\lvert H^{[0]}(w)\rvert\le C_\varepsilon e^{(2/\pi)\lvert\operatorname{Im}w\rvert+\varepsilon\lvert w\rvert}$ stands and is what makes Lemma 7's strips work.

**R5 (cosmetic).** Lemma 4: state the $m=0$ convention ($\lvert H_0\rvert=2\le C_H$). Theorem 2 Tonelli: "$\lvert\eta\rvert^j/j!^2$" should be "$e^{\lvert\eta\rvert}\lvert\eta\rvert^j/j!$" (harmless). Gevrey paragraph: the $N!$ comes from the Laplace integral, not from $R_N$'s Cauchy bound. Theorem 3: "$B-H$'s nearest singularities are at distance $2\rho$" → at $i\pi$ (distance $\pi/2$); conclusion unaffected.

## 3. Evidence table (all printed by referee's code; digit counts are residuals, not claims)

| Check | File | Result |
|---|---|---|
| $E^*=5/24$ certification | step_03, step_08, verify[1] | piece A $c_1+(c_2+T)u\le-1.7885\mathrm e{-3}<0$; piece B $\max\psi=0.0433847019$; piece C $0.0399729017$; all $<(5/24)^2$ |
| two-term law vs exact $H_n$ | step_04, verify[2] | $\lvert H_n-\mathrm{law2}\rvert\le4.58\rho^{-n}n^{-7/2}$ measured; $C_H=5.8739238$ at $n=3$; $\min\lvert H_n/\mathrm{lawH}\rvert=0.87454216$ at $n=10$ |
| claimed $4.45$-bound truth | step_11, verify[2] | worst $\lvert R_n\rvert/(4.45u+135(2.9\rho)^{-n})=0.9864$ ($n=18$), $n\le220$ — TRUE |
| repaired $7.06$-bound | step_07, 09, 10, verify[3] | certified sup $7.05229299615373$ ($n=11$) on $10..500$; Stirling tail $\le5.19$ ($n\ge500$) |
| false $1.63$-intermediate | step_33 | $\Gamma(7.5)/\Gamma(11)=1.6306839\cdot10^{-3.5}>1.63\cdot10^{-3.5}$ |
| **D** bound vs true tail | step_12, verify[4] | $n=20/40/80/120$: $\lvert\mathrm{tail}\rvert=2.84\mathrm e5/1.78\mathrm e{24}/3.45\mathrm e{74}/1.50\mathrm e{134}\le$ bound $9.4\mathrm e8/1.8\mathrm e{28}/1.0\mathrm e{79}/8.3\mathrm e{138}$; relative to main $9.2\mathrm e{-8}, 3.9\mathrm e{-14}, 1.7\mathrm e{-26}, 1.0\mathrm e{-38}$ |
| **E** per-class Richardson | step_13, verify[5] | const $0.9999977$–$0.9999978$ (all 4 classes); $c_1=-0.80202$ (pred $(3-3\pi)/8=-0.80310$) even, $+1.55414$ (pred $+1.55310$) odd; residual $\times n^2\to-4.05$ (even), $+0.43$ (odd) |
| Lemma 2 (B) | step_22–24, verify[7] | $\lvert\sum H_{n-j}M_{n,j}-f_n\rvert/\lvert f_n\rvert\le1.84\mathrm e{-59}$, $n\le8$; $M_{1,1}=(\ln2-1)/2$ exact; $M$-bound worst ratio $1.0$ (tight at $j=1$) |
| Claim A spot | step_28, verify[8] | Laplace form vs perimeter at $\eta=0.37$: $\lvert\mathrm{diff}\rvert=7.8\mathrm e{-62}$ |
| Lemma 5/6 identities | step_14, 15, 19 | $J_0$ form vs defining series $\sim1\mathrm e{-8}$ ($\mathbb R^+$); Bessel-$I$ form $\sim1\mathrm e{-7}$ ($\mathbb R^-$) |
| Lemma 7 $\mathfrak B_j$ decomposition | step_19–21 | disk $z=\pm1$: $2.7\mathrm e{-12}, 2.1\mathrm e{-10}$; off disk $z=+4$: $4.0\mathrm e{-10}$, $z=-4$: $4\mathrm e{-9}$ (vs independent conformal map) |
| Lemma 8 identity | step_25–26 | $\int_0^\infty e^{-t}H^{[j]}(wt)dt=w^jH(w)$ at $w=\pm0.5,\pm2,-5$, $j=0,1,2$: diffs $1\mathrm e{-10}$–$2\mathrm e{-5}$, on/off disk |
| **Theorem 2**: $S(\eta)=I(\eta)$, own Borel–Laplace | step_16–17, 32, verify[6] | conformal map $\zeta=2\pi w/(1-w^2)$ (exact rational map), H-split; vs perimeter: $\eta=0.5$: $3.5\mathrm e{-24}$; $\eta=2$: $8.7\mathrm e{-10}$; $\eta=-0.5$: $4.6\mathrm e{-23}$; $\eta=-0.9$: $4.8\mathrm e{-12}$ |
| G: radius of $B-H$; $T2b$ decomposition | step_16, 27 | $\lvert r_{180}\rvert^{1/180}=0.30484\to1/\pi=0.31831$; $R_2=\sum_{j\ge1}\mathfrak B_j$ vs map-side: $6.6\mathrm e{-25}$ ($z=2$), $5.0\mathrm e{-22}$ ($z=2.5$), $1.0\mathrm e{-7}$ ($z=3$) |
| H jump asymptotics | step_29, verify[9] | $\lvert J\rvert/\lvert\mathrm{pred}\rvert=0.958,0.981,0.989,0.995,0.997,0.9988$ at $\varepsilon=0.3..0.05$ ($\to1$); phase defect $\approx-3\varepsilon/4$; $J/\mathrm{pred}-1=1.2\mathrm e{-6}+1.0\mathrm e{-4}i$ at $\varepsilon=0.05$; $O(\varepsilon^2)$ coefficient extracted $\to-0.78125=-25/32$ |

Also: sign pattern $\operatorname{sign}(f_n)=(-,+,+,-)$ for $6\le n\le59$ measured (deviations only $n=0,1,4,5$), provable for all $n\ge10$ from the two-term bounds; Lemma 3.4's nonvanishing-orbit argument spot-checked for three rhombi (step_30).

## 4. Paper-sentence answer

1. **"(We have a Watson–Darboux argument, but one uniformity step is not written out, so we state the law as numerical.)"** → The uniformity step is now written out and proved (Lemma 4, explicit constants, verified). **The growth law — square two-term law with $O(n^{-2})$, rhombus law, and thm:zero's radius zero — may be called proved**, with a † sketch in the paper and the full proof in the notes (after repairs R1–R3, which change no displayed formula of GRW1 except the proof-internal remainder constant $4.45\to7.06$).
2. **"…reproduces direct quadrature of $I(\eta)$ at $\eta=0.1,0.2,0.3$ and $-0.2$ to 44, 31, 25 and 30 digits."** → The first three reproduce (44.6/31.3/25.4 with $[90/90]$); **the "30" at $\eta=-0.2$ is not reproduced** — plain Padé–Borel gives 19.3/21.0/23.8 digits at $m=60/75/90$ (planner; GRW1 concurs). Change 30 → 24, or attribute the figure to the H-split variant (42.7 digits). The Borel summability justifying the resummation is now proved (Theorem 2, both real directions, every $\eta>-1$, after repairs R4–R5).
3. The Richardson sentence ($6\cdot10^{-6}$, $n\le64$) is consistent with the referee's $2.3\mathrm e{-6}$ at $n\approx152$–$176$.

Deliverables in `new_results/GRW1R_referee_diamond/`: `result.md` (full report with line-by-line ledger), `verify_grw1r.py` (52 s, reproduces the whole evidence table), `progress.log`, `step_01.py … step_36.py`.