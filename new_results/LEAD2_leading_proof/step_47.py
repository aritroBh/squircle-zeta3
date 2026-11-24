import os, time
os.chdir('<project>/research/08_theory/LEAD2_leading_proof')
full = open('proof_lead.md').read()
old1 = """**Lemma 2.1 (exact algebra).** With `u = (v-1/v)/2`, `v = u + S` (positive branch, `v > 0` for
`u in (-1,1)`... `v` runs from `b` to `a` as `u` runs from `-1` to `1`):"""
new1 = """**Lemma 2.1 (exact algebra).** With `u = (v-1/v)/2`, `v = u + S` (the branch with `v > 0`;
as `u` runs from `-1` to `1`, `v` runs from `b` to `a`):"""
assert old1 in full; full = full.replace(old1,new1)
old2 = """**Proof.** Rational-function identities in `Q(sqrt2)(u)[S]/(S^2 = 1+u^2)`, reduced symbolically
(sympy `radsimp` gives 0 for each difference; `verify_lead2.py` part B). For (c): `(1-u^2) m_e` even,
and `(1-u^2) m_e(u) = 1 - sqrt2/S -> 0` at `u = ±1` with rate `(u∓1)(±1)/4·...`, giving `-1/4`;
(2) is `(1+u - u + u sqrt2/S)/(1-u^2) = m_odd`.  ∎"""
new2 = """**Proof.** Rational-function identities in `Q(sqrt2)(u)[S]/(S^2 = 1+u^2)`, reduced symbolically
(sympy `radsimp` gives 0 for each difference; `verify_lead2.py` part B). For (c): `(1-u^2) m_e` is
even in `u`, and `1 - sqrt2/S = (S^2-2)/(S(S+sqrt2)) = (u^2-1)/(S(S+sqrt2))`, hence
`m_e(u) = -1/(S(S+sqrt2)) -> -1/4` at `u = ±1` (in particular `m_e` is continuous on the closed
interval); (2) is `(1+u) - u(1 - sqrt2 u·... )`: directly, `1/(1-u) - u(1-sqrt2/S)/(1-u^2)
= (1 + u sqrt2/S)/(1-u^2) = m_odd`.  ∎"""
assert old2 in full; full = full.replace(old2,new2)
old3 = """with `M(i,j; m) := int_{-1}^1 L1^i L2^j m(u) du` (convergent for every occurring pair: near `u = 1`,
`m_odd ~ 1/(1-u)` but `L1 ~ -(1-u)·...` supplies convergence since `i >= 1` throughout; near `u = -1`
all measures are finite except `m_e, m_oddlast` which are still `O(1)`, and `L1^i` is integrable):"""
new3 = """with `M(i,j; m) := int_{-1}^1 L1^i L2^j m(u) du` (convergent for every occurring pair: at `u = -1`
all four measures are finite and `L1^i` has at worst an integrable log singularity; at `u = 1`,
`m_e, m_diag` are finite by (c), while `m_odd ~ 1/(1-u)`, `m_oddlast ~ 1/(2(1-u))` — but then
`L1 = log((1+u)/2) = O(1-u)`, and every occurring key has `i = d-j >= 1`, so
`L1^i L2^j m_odd = O((1-u)^{i-1} L2^j)` is integrable):"""
assert old3 in full; full = full.replace(old3,new3)
open('proof_lead.md','w').write(full)
p2 = open('proof_part2.md').read()
open('proof_lead.md','w').write(full)
print('patched section 2')
with open('progress.log','a') as f:
    f.write(f"{time.strftime('%H:%M:%S')} step: cleaned section 2 artifacts (Lemma 2.1 proof, convergence note)\n")
