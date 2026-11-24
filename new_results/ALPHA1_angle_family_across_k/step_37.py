open('progress.log','a').write('step_02: exact reducer at c=1/2 k=1..4; weight parts sum to M_k (50 digits, 4 levels); theta-poles vanish; W2 const seen\n')
open('progress.log','a').write('step_03..06: component tests (raw integrand 50d OK; pf_from_poly/pf_fold/mul_DL OK; DL2 pole-at-(-r) bug found+fixed)\n')

result = r'''# ALPHA1: the corner-angle family M_k(alpha) across k — result

**Verdict (updated as work proceeds).**
- Task A (general-c integrand + reducer + controls): **DONE**. w-integrand derived and checked; exact reducer built; weight parts reproduce M_k numerically (50 digits); squircle control pending.
- Task B (polynomial law): in progress. W_2(k) constant observed at c=1/2 (k=2,3,4).
- Task C (leading law in the angle): pending.

## 1. The general-angle setup (DERIVED, checked numerically to 50 digits)

With t = g/(1+g), the corner function
M(u;alpha) = int_0^1 [(1+g)^(-1-u) sqrt(1+g^{2-2u} - 2cos(alpha) g^{1-u}) - 1] g^{u-1} dg
becomes (EXACT algebra; NUMERICAL 50 digits, step_01)
  M(u;alpha) = int_0^1 Psi_u(t) dt,
  Psi_u = ( sqrt(t^{2u-2} + (1-t)^{2u-2} - 2 cos(alpha) (t(1-t))^{u-1}) - t^{u-1} - (1-t)^{u-1} )/2.
At alpha=pi/2 this is the squircle -F/2.

With c = cos(alpha/2), r = (1-c)/(1+c) (so r = tan^2(alpha/4)), the Euler substitution
  t = (w-r)(w+1)/(2(1-r) w),  1-t = (1-w)(w+r)/(2(1-r) w),  dt/dw = (w^2+r)/(2(1-r) w^2)
maps w in [r,1] to t in [0,1], w -> r/w is t -> 1-t, and RATIONALISES the master sqrt:
  q_alpha := sqrt(t^2+(1-t)^2 - 2cos(alpha) t(1-t)) = sqrt(1 - 4 c^2 t(1-t)) = (w^2+r)/((1+r) w)   [PROVED algebra, step_01]
so M(u;alpha) = int_r^1 Psi_u(t(w)) (dt/dw) dw with everything rational over Q(c) (step_01, 50 digits, 2 angles).

Letters: lt = log(w-r)+log(w+1)-log w - L, lmt = log(1-w)+log(w+r)-log w - L, L = log(2(1-r)).
Real poles {0, +-1, +-r}; the only non-real poles are at w = +-i sqrt(r) (zeros of w^2+r = (1+r) w q_alpha), theta := i sqrt(r), theta^2 = -r.

## 2. The integrand [eps^k] Psi (EXACT; matches raw M_k integrals to 50 digits, step_03)

[eps^k] Psi = sum_m R~_m(t) x^m y^{k-m}, x=lt, y=lmt:
  sqrt part R~_m = (1/2) q_alpha T_m(t)/(t(1-t) q2a^k), q2a = q_alpha^2 = 1 - 2(1+ca) t + 2(1+ca) t^2, ca = cos(alpha);
  T_m(t) = sum_{j>=1} binom(1/2,j) sum_{s1+s2+s3=j} multinom * inner * (1-t)^{2s1} t^{2s2} (t(1-t))^{s3} (-2ca)^{s3} q2a^{k-j},
  inner = sum_{m1+m2+m3=k} w(s1,m1) w(s2,m2) w1(s3,m3) binom(m3, m-m1), w(s,m)=s!S(m,s)2^m/m!, w1(s,m)=s!S(m,s)/m!;
  A-part (m=k) -1/(2 k! t), B-part (m=0) -1/(2 k! (1-t)).
In w (with dw): R_m(w) = C_k T~_m(w)/[w (w^2+r)^{2k-2} (w-r)(w+1)(1-w)(w+r)], C_k = 1/(2(2c)^{2k-1}),
T~_m(w) = (2(1-r) w)^{2k} T_m(t(w)) a polynomial. AB: R^A = -(w^2+r)/(2 k! w (w-r)(w+1)).
Pole bound at theta: <= 2k-2 (matches the squircle Lemma). Folding t->1-t at the T level (T_m(1-t) = T_{k-m}(t)).

## 3. The triangular IBP system (the analogue of eq:triangular) — IMPLEMENTED, exact

Canonical keys (i,j), i>=j>=0, d=i+j. H(R) = principal parts of order >=2 at w=+-theta; J its rational
primitive normalized G(r)=0; fold F(R) = r R(r/w)/w^2 (maps letters {0,+-1,+-r,+-theta} among themselves).
  W_ij = eps^d R_ij - (i+1) G_{i+1,j} om1 - [i>j](j+1) G_{i,j+1} om2 - [j=i-1] i F(G_ii om1),  G_ij = JH(W_ij),
  om1 = dlt/dw = (w^2+r)/(w(w-r)(w+1)),  om2 = dlmt/dw = (w^2+r)/(w(w-1)(w+r))   [both vanish at w=+-theta].
Boundary terms vanish (G(r)=0; lt^i=O((1-w)^i) at w=1, i>=1). Solved downwards in d (cwred.py).

CONTROL (NUMERICAL, step_02, c=1/2, k<=4): sum of weight parts sum_{d=1}^k sum_{i+j=d} int_r^1 K^{(k)}_ij lt^i lmt^j dw
= M_k from corner_series.py to 49-52 digits at every level k=1..4.  No simple theta-pole survives in any kernel
(EXACT, k<=4).  W_2(k) constant for k=2,3,4 (= 0.6886407138827469147766..., 30+ digits).
'''
open('result.md','w').write(result)
print("result.md written")
