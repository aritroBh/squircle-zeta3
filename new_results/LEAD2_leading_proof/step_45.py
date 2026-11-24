import os, time
os.chdir('<project>/research/08_theory/LEAD2_leading_proof')
p2 = open('proof_part2.md').read().rstrip() + "\n\n"
thm = """**Theorem 6.2 (the law).** Given the D-table of Section 1 for all `d`, for every `w >= 2`:

    c(w) = eps_{w-2} A_{w-1}/(2 (w-2)!) * e_w,   e_w = [y^w] E(y).

**Proof.** `c(d+1) = delta_d K_d/(d-1)!` (definition of `K_d`) `= delta_d (d!/2) e_{d+1}/(d-1)!`
(Theorems 4.2 and 5.5) `= delta_d (d/2) e_{d+1}`. With `delta_d = eps_{d-1} A_d/d!` (Lemma 6.1) the
prefactor is `eps_{d-1} A_d/(2 (d-1)!)`, and `w = d+1` gives the statement. ∎

**Corollary 6.3 (generating functions).**
  (a) `F(u) := sum_{w>=2} c(w) u^w = (1/2)[ u^2 (1-sinh u)/cosh^2 u  ⊙  E(u) ]`
      (Hadamard = coefficientwise product), since `eps_{w-2} A_{w-1}/(w-2)! = [u^{w-2}]
      (u^2 (1-sinh u)/cosh^2 u)` by Lemma 6.1.
  (b) The reduced GF is a pure Gamma product:
      `sum_{w>=2} [2(w-2)!/(eps_{w-2} A_{w-1})] c(w) u^w = E(u) - 1 - u log 2`.
  (c) `E(y) = 2^y Gamma(1-y) Gamma(1+y/2)/Gamma(1-y/2)`: Legendre duplication
      `Gamma(z)Gamma(z+1/2) = 2^{1-2z} sqrt(pi) Gamma(2z)` at `z = 1/2 - y/2` gives
      `Gamma(1/2-y/2) = 2^y sqrt(pi) Gamma(1-y)/Gamma(1-y/2)`, substitute in the definition of `E`.
      Hence with `log Gamma(1+t) = -gamma t + sum_{n>=2} (-1)^n zeta(n) t^n/n`,
      `log E = y log 2 + sum_{n>=2} s_n y^n/n`, `s_n = zeta(n)` (n even), `s_n = (1-2^{1-n}) zeta(n)`
      (n odd).  [All checked symbolically/numerically, part G.]
"""
p2 = p2 + thm
open('proof_part2.md','w').write(p2)
p1 = open('proof_part1.md').read()
# rebuild app from existing proof_lead.md section 7
oldfull = open('proof_lead.md').read()
app = oldfull[oldfull.find("\n## 7. Verification guide"):]
full = p1 + p2 + app
open('proof_lead.md','w').write(full)
print('sections:', [s for s in ['## 1.','## 2.','## 3.','## 4.','## 5.','## 6.','## 7.'] if s in full])
print('has Theorem 6.2:', 'Theorem 6.2' in full, ' Corollary 6.3:', 'Corollary 6.3' in full)
print(len(full))
with open('progress.log','a') as f:
    f.write(f"{time.strftime('%H:%M:%S')} step: restored Theorem 6.2 + Corollary 6.3 (lost in earlier edit); reassembled proof_lead.md\n")
