PRV1 complete — the depth-2 reductions of the regular octagon and decagon are now **proved** (exact symbolic), and the general mechanism is isolated. Work directory: `research/08_theory/PRV1_depth2_proofs/` (full write-up `result.md`, verifier `verify_full.py` → `ALL CHECKS PASSED`, log `verify_full.out`, trail `step_01.py`–`step_66.py`, `progress.log`).

## Theorem O (octagon; CX1's verifier made fully exact + OCT6 matched exactly)
With $\operatorname{Li}_{2,1}(x,y)=\sum_{n>m\ge1}x^ny^m/(n^2m)$, $\Lambda_{\rm oct}=\operatorname{Li}_{2,1}(1-\sqrt2,-\sqrt2/2)$ and $\mathrm{swap}_{\rm oct}=\operatorname{Li}_{2,1}(1-\sqrt2/2,-\sqrt2)$ equal:
- $\Lambda_{\rm oct}=\frac34\zeta(3)+\operatorname{Li}_3(k)+2\operatorname{Li}_3(x_1)-\operatorname{Li}_3(d)+\frac{l_2}{2}\operatorname{Li}_2(k)+\frac{l_2}{4}(\operatorname{Li}_2(x_2)-\operatorname{Li}_2(h))$
- $\mathrm{swap}_{\rm oct}=\operatorname{Li}_3(d)+2\operatorname{Li}_3(x_2)-\operatorname{Li}_3(k)-4\operatorname{Li}_3(h)-4\operatorname{Li}_3(-h)-\frac{l_2}{2}(\operatorname{Li}_2(k)+\operatorname{Li}_2(x_2))+\frac{l_2^2}{8}(\frac{l_2}{2}-l_u)$

($k=\sqrt2-1$, $x_1=1-\sqrt2$, $x_2=1-\sqrt2/2$, $d=2-\sqrt2$, $h=\sqrt2/2$, $l_2=\log2$, $l_u=\log(1+\sqrt2)$.) These equal OCT6's stored `reductions.pkl` **exactly** (symbolic difference 0, a rigid comparison).

## Theorem D (decagon; previously 830-digit LLL conjectures, now proved)
$\Lambda_{\rm dec}=\operatorname{Li}_{2,1}(-1/\varphi,-1/\varphi)=\frac25\zeta(3)+\frac45 l_\varphi z_2-l_\varphi^3+3\operatorname{Li}_3(-1/\varphi)+\operatorname{Li}_3(1/\varphi)$ and $\mathrm{swap}_{\rm dec}=\operatorname{Li}_{2,1}(1/\varphi^2,-\varphi)=\frac35\zeta(3)-\frac75 l_\varphi z_2+2l_\varphi^3+\operatorname{Li}_3(-1/\varphi)+\operatorname{Li}_3(1/\varphi)$ — exact match with DEC1 modulo the proved golden sum relation $\operatorname{Li}_3(1/\varphi)+\operatorname{Li}_3(-1/\varphi)=(\zeta(3)-l_\varphi z_2)/5+l_\varphi^3/6$.

## Proof method (both, self-contained)
(a) **Foundations**: $\partial_x\operatorname{Li}_{2,1}=\operatorname{Li}_{1,1}/x$, $\partial_y\operatorname{Li}_{2,1}=(\operatorname{Li}_2(x)-\operatorname{Li}_2(xy)/y)/(1-y)$ (termwise on $|x|,|y|<1$; checked as exact finite-sum identities), and $\operatorname{Li}_{1,1}(x,y)=\log(1-x)\log(1-y)+\operatorname{Li}_2(\frac{-y}{1-y})-\operatorname{Li}_2(\frac{-y(1-x)}{1-y})$ for $-1<x,y<1$ (proved by $\partial_x$ + pin at $x=0$, exact). (b) One-parameter $C^1$ families on $0<t<1$: octagon $G(t)=\operatorname{Li}_{2,1}(-t,-\frac{1-t}{2t})$, $H(t)=\operatorname{Li}_{2,1}(\frac{1-t}{2},-\frac{2t}{1-t})$; decagon $G_{\rm dec}(t)=\operatorname{Li}_{2,1}(-t,-\frac{1-t}{t})$, $H_{\rm dec}(t)=\operatorname{Li}_{2,1}(1-t,-\frac{t}{1-t})$; chain-rule derivatives assembled exactly from the foundations (20 rational curve identities, residuals 0). (c) Explicit single-polylog closed forms $G_{\rm cl},H_{\rm cl},G_{{\rm dec},{\rm cl}},H_{{\rm dec},{\rm cl}}$ (in result.md) whose derivatives equal those of $G,H,G_{\rm dec},H_{\rm dec}$ **exactly**, modulo classical dilog identities E (Euler), L1/L2 (Landen), DUP2/DUP3 (duplication), OC (odd complement) — each proved by differentiation (exact residual 0) + exact endpoint pin. (d) Constants pinned exactly: $G,G_{\rm dec}\to0$ at $t\to1^-$ with $G_{\rm cl}(1)=G_{{\rm dec},{\rm cl}}(1)=0$ (uses $\operatorname{Li}_3(-1)=-3\zeta(3)/4$); $H,H_{\rm dec}\to0$ at $t\to0^+$ with $H_{\rm cl}(0)=H_{{\rm dec},{\rm cl}}(0)=0$. (e) Specialisation at $t=\sqrt2-1$ (over $\mathbb Q(\sqrt2)$) and $t=1/\varphi$ (over $\mathbb Q(\sqrt5)$, golden lattice incl. Lewin's $\operatorname{Li}_3(1/\varphi^2)=\frac45(\zeta(3)-l_\varphi z_2)+\frac23l_\varphi^3$), all argument identities exact.

## Lemma L (general result)
For a Landen pair $x_2=x_1/(x_1-1)$ (two-ratio edges of regular $2m$-gons), set $c=x_2/(1+x_1)$. The line $p=c(1+x)$ through $(x_1,x_2)$ gives $G_c,H_c$ with $G_c(t)\to0$ ($t\to1^-$), $H_c(t)\to0$ ($t\to0^+$); for **every** $c>0$ the derivatives reduce to depth 1 (Lemma 2 applies: the Möbius map $T_c(t)=c(1-t)/(c+(1-c)t)$ is an involution, checked exactly); the primitive closes in single polylogs for $c=\frac12$ (octagon, $m=4$) and $c=1$ (decagon, $m=5$) — the two proved cases. Square ($m=2$) has no depth-2 constants; the $c=\frac23$ member (hexagon-type pair) shows **no** classical closure (PSLQ artifact test below), and the hexagon's actual depth-2 content is Nielsen $S_{1,2}$-type, already reduced by K1e/K1f.

## Evidence table (all numbers printed by code; exact parts residual ≡ 0 symbolically)
| Check | dps | residual |
|---|---|---|
| verify_full.py: Parts 0,1,CR,2,2b,3,3b,3c (all exact sympy) | — | **0 exactly; ALL CHECKS PASSED** |
| OCT6 `reductions.pkl` vs §1.5 forms (both) | exact | 0 |
| DEC1 forms vs specialisation (both) | exact | 0 |
| Classical ids E,L1,L2,DUP2,OC,DUP3,Landen-3 (3 pts) | 120 | ≤3.9e-121 |
| Lemma 2 vs series (3 pts) | 120 | ≤5.9e-122 |
| $G{=}G_{\rm cl}$, $H{=}H_{\rm cl}$; $G_{\rm dec}{=}$, $H_{\rm dec}{=}$ (5 pts each) | 120 | ≤1.9e-121 |
| $\Lambda_{\rm oct}=-0.0275883613456533384687572087639\ldots$ vs $G_{\rm cl}(k)$ | 120 | 1.7e-121 |
| $\mathrm{swap}_{\rm oct}=-0.0321696900151107721899470611759\ldots$ vs $H_{\rm cl}(k)$ | 120 | 1.6e-121 |
| $\Lambda_{\rm dec}=-0.05126163058647632738126326036571123343019\ldots$ vs DEC1 | 170 | 4.5e-172 |
| $\mathrm{swap}_{\rm dec}=-0.06341782168216461215085970440865279340954\ldots$ vs DEC1 | 170 | 1.2e-171 |
| Golden lattice (Lewin, sum relation, 3 dilogs) | 170 | ≤5.2e-172 |
| $c=\frac23$ scope test: PSLQ candidate (coeffs ~5e7, residual 2.2e-110) | 150→250 | frozen at 2.2e-110 ⇒ spurious; no closure |
| CX1's original verify.py re-run | — | ALL CHECKS PASSED (its H pin was numeric, -7.1e-102; ours is exact) |

Key files: `verify_full.py` (self-contained exact verifier + numeric evidence), `verify_full.out`, `result.md` (full proofs §1 octagon, §2 decagon, §3 general lemma, §4 verification), `progress.log`, `step_01.py`–`step_66.py`.