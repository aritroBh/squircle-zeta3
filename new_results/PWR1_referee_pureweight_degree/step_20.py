log = open('progress.log','w')
log.write("PWR1 referee: thm:pw (pure weight) + prop:degree (Fourier coefficient degree)\n")
log.write("step_01-14: exploration: read paper main.tex (sec:pure, sec:fourier), D16/result.md, H1/result.md, K3/result.md, K3 step_12/13/16, e_rm.pkl\n")
log.close()
result = open('result.md','w')
result.write("""# PWR1 referee report: thm:pw (pure weight) and prop:degree (degree of Fourier coefficients)

**Referee:** independent re-derivation + exact-symbolic and high-precision numerical tests.
**Objects.** (A) Theorem thm:pw: log(P(p)/(4 sqrt 2)) = sum_{k>=2} w_k tau^k, w_k in Q[zeta(2),..,zeta(k)]
homogeneous of weight k, tau = 1-1/p; same for rhombi; with rational prefactor Rat_n(tau) for n-dim surface area.
Written proof under review: research/08_theory/D16_pure_weight/result.md (Proof (mechanism) M1-M3,
single-appearance lemma, d16_engine.py) + H1 result.md sec.3 (T6).
(B) Proposition prop:degree: [tau^n] (-1)^k (pi/sqrt2) a_k(h_p) = Theta_n(k^2), deg exactly floor(n/2)-1,
alternating-MZV coefficients of weight n; leading coeff for n=2m: (-1)^{m-1}(2^{2m-1}-1) zeta(2m).
Written proof under review: research/08_theory/K3_fourier_all_orders/result.md sec.1.1-1.3.

## Referee's re-derivation (summary; details below)

### A. thm:pw
[Re-derived: H1 T1 representation == D16 B_n representation via Lukacs/size-bias shift
 Gamma(1+a) -> Gamma(eps), 1+a-tau/n = eps (checked exactly). M1 moment formula re-derived:
 E[prod (e^{2 tau m_i}-1)] = sum_S (-1)^{j-s} prod_v Gamma(1+a+2 tau(nu_v-s/n))/Gamma(1+a)^n -- OK.
 gamma cancels in each log-moment (sum_v (nu_v - s/n) = 0) and in the flat-face factor
 (n(-gamma a) = +gamma(n-1) tau vs -log Gamma(1-(n-1)tau) = -gamma(n-1)tau) -- OK.
 Pure weight: every log Gamma argument is tau x rational, so zeta(k) enters at exactly tau^k;
 exp/products/finite sums preserve grading; B_n|_{Z=0}=1 via sum_s binom(j,s)(-1)^{j-s} n^{j-s} n^s = n^j(1-1)^j = 0 -- OK.
 Rat_n(tau) = (1-tau)^{n-1}/prod_{i=1}^{n-1}(1-(n-1)tau/i) = prod_{i=1}^{n-2} i/(i-(n-1-i)d) -- re-derived OK.
 Lambda_n = -sum zeta(k)(n-1)^k(1-n^{1-k}) tau^k/k -- re-derived OK.]

### B. prop:degree
[Re-derived: a_k - A_k = (-1)^k B(tau,2k), B even entire in kappa; layer identity
 h_p - h_1 = (1+e^{-2y})^{-1/2}(e^{tau L(y/tau)}-1), y=-log tan theta, L(s)=log(1+e^{-s});
 Watson substitution y = tau s gives D(tau;m) = (8/pi) tau int H(tau s) gd(tau s)^{2m} sum_j tau^j L^j/j!;
 D~_{n,m} = 2 sum_{j=1}^{n-1-2m} e_{m,n-1-j} I_{n-1-j,j}/j!, e rational, I alternating MZV weight n;
 degree <= floor(n/2)-1 since j>=1; m-sum at fixed tau^n is FINITE (m <= n/2-1), so no uniformity-in-k issue;
 tops: D~_{2n,n-1} = 2(2n-2)! eta(2n) (only j=1; e_{n-1,2n-2}=1; I_{2n-2,1}=(2n-2)! eta(2n)) != 0;
 odd: D~_{2n+1,n-1} = (2n-1)! eta(2n+1) + 2(2n-2)![eta(2n+1)-A_{2n}] with classical A_{2n}
 = (2n-2)![zeta(2n+1)+2 sum_{j=1}^{n-1} eta(2j) zeta(2n+1-2j)] > 0. Leading-coefficient relation
 2(-4)^{m-1}(1-2^{1-2m}) = (-1)^{m-1}(2^{2m-1}-1) checked algebraically.]

## Verdict tables and numerical evidence: TO BE FILLED
""")
result.close()
print("wrote progress.log, result.md skeleton")
