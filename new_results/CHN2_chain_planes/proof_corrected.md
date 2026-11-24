# CHN2: corrected elementary proof (referee workflow wf_add8df88-255, judge text; verdict: CORRECT WITH FIXES. The mathematics is correct as written. I...)

Referee verdict: CORRECT WITH FIXES. The mathematics is correct as written. I re-derived Steps 1 to 3 and the final algebra by hand and found no error. The rigor and blind lenses agree, and the numeric attack found no counterexample. I also ran an independent exact rational check in pure Python. The claim holds with zero mismatches at all 128 integer plane points (a,b in 1..8 on each of Pi1 and Pi2). It reproduces F23(2,2,1,1)=1/15 and F23(2,1,2,1)=4/45. It fails off both planes, with F23 = 5/84 against a formula value of 11/168 at (2,1,3,1), and 2/315 against 1/140 at (1,2,3,4). The one substantive omission is in the writing, not the mathematics. Step 2 applies Fubini and the fundamental theorem of calculus without checking integrability or the corners, which matters because the integrands blow up at the edges and corners when some u_i < 1. The slice argument in the corrected proof closes this, and it needs no condition beyond u_i > 0. The other fixes are notational and scope statements. The plane hypothesis is used only in Step 3. Step 2, F23 - F32 = G0, holds for all u > 0. The proof is self-contained and uses none of the cited papers.

Required fixes (all incorporated below):
- Step 2 (the only substantive fix): after the divergence identity, add the integrability and slice justification. A, X and C below are products of one-variable factors s^a(1-s)^b with a,b > -1, so they lie in L^1 of the unit square for all u > 0. For each fixed t in (0,1), the map s -> P(s,t) is continuous on [0,1-t], is C^1 on (0,1-t], vanishes at s=0 because u2 > 0, and has an integrable derivative, so its s-integral equals P(1-t,t) exactly. The same holds for Q in t with s fixed, using u3 > 0. State that no corner terms arise, because no limit along shrinking arcs is taken and the corners are a null set. A closed-triangle divergence theorem would not apply when u3 < 1 or u1+u3 < 1. State that the two hypotenuse integrals converge because u1+u3 > 0.
- Notation: the letter B currently means the body, the unit-square integral, the monomial B' of Step 2, and the Beta function. Rename the unit-square integral to V and the monomial B' to C, and keep Beta only for the Beta function. Also name X := s^u2 (1-s)^(u1-1) t^(u3-1) (1-t)^u4, the cross term that cancels. Write Gamma consistently.
- Scope statement: say that the claim is proved for real u_i > 0 on Pi1 or Pi2. Say that Step 2 holds on all of (0,inf)^4 and that only Step 3 uses the plane hypothesis. Do not state or imply that Pi1 and Pi2 are the only loci where the formula holds. That is not proved, and numerics show it fails off the planes. If a complex extension is wanted, state it only for Re u_i > 0 with the plane relations read as complex equalities (see analytic_continuation_domain).
- Step 3: say explicitly that (s,t) -> (1-t,1-s) is a linear involution of the unit square with |det| = 1, mapping {s+t >= 1} onto T because s'+t' = 2-(s+t). Also say that the diagonal s+t=1 is a null set, so B = F23 + I_up in Step 1 is exact.
- Preamble: keep F23 = Vol(B)/16 as a remark with the one-line derivation (see the corrected proof). It lies outside the claim as stated. Only the rigor lens verified it, by hand and by Monte Carlo. The blind lens did not check it.
- Citations: use the exact in-text wording in the citations list. Always say 'arXiv v3 numbering' (1205.1516), 'v2 numbering' (1304.7267 and 1304.7304), or 'v2' (1908.10830) when giving an equation number. Describe the P Q exp form as an observed pattern tested to weight 16, not a theorem. Quote the weight-8 basis whenever the coefficient 1/5 of zeta_{3,5}[M5,M3] is quoted. Do not write 'F23 = F^(23) of Schlotterer-Stieberger' unless the identification is separately established (see residual_doubts).

## Corrected proof

CLAIM. Let u1,u2,u3,u4 > 0 and phi(s,t) = u2 u3 s^(u2-1) t^(u3-1) (1-s)^u1 (1-t)^u4. Let T = {s,t >= 0, s+t <= 1} and F23 = integral of phi over T. If (u1=u2 and u3=u4) [plane Pi1] or (u1=u3 and u2=u4) [plane Pi2], then
 F23 = (1/2)[ G0' + G0 ], where G0' = G(1+u1)G(1+u2)G(1+u3)G(1+u4) / (G(1+u1+u2)G(1+u3+u4)) and G0 = G(1+u1+u3)G(1+u2+u4)/G(1+u1+u2+u3+u4), with G = Gamma.

Remark (not part of the claim). With B = {|x1|^p1+|x2|^p2 <= 1, |x2|^p2+|x3|^p3 <= 1, |x3|^p3+|x4|^p4 <= 1} and u_i = 1/p_i, take the orthant x_i >= 0 and put x2 = s^u2, x3 = t^u3. Then dx2 dx3 = u2 u3 s^(u2-1) t^(u3-1) ds dt. The fibre of x1 has length (1-s)^u1, the fibre of x4 has length (1-t)^u4, and the middle constraint reads s+t <= 1. There are 16 orthants, so Vol(B) = 16 F23.

NOTATION. Beta(x,y) = G(x)G(y)/G(x+y). Define the one-variable-product functions
 A = s^(u2-1) t^(u3-1) (1-s)^u1 (1-t)^u4, so phi = u2 u3 A;
 X = s^u2 (1-s)^(u1-1) t^(u3-1) (1-t)^u4;
 C = s^u2 (1-s)^(u1-1) t^u3 (1-t)^(u4-1).
Define V = integral of phi over [0,1]^2, I_up = integral of phi over {s+t >= 1} within [0,1]^2, and F32 = u1 u4 * integral over T of C.

Integrability. A, X and C are each a product of a function of s and a function of t, and each factor has the form s^a (1-s)^b. The exponents that occur are u2-1, u1, u2, u1-1, u3-1, u4, u3, u4-1, all > -1 because every u_i > 0. Hence A, X and C are integrable on [0,1]^2, and so on T. All integrands are nonnegative, so Tonelli applies freely. In particular F23 and F32 are finite.

Step 1. The diagonal s+t=1 is a null set, so V = F23 + I_up. The integrand factorises, so
 V = [u2 Beta(u2,1+u1)] [u3 Beta(u3,1+u4)] = G(1+u1)G(1+u2)/G(1+u1+u2) * G(1+u3)G(1+u4)/G(1+u3+u4) = G0',
using x Beta(x,1+y) = G(1+x)G(1+y)/G(1+x+y).

Step 2 (all u > 0, no plane needed): F23 - F32 = G0.
Put P(s,t) = u3 s^u2 (1-s)^u1 t^(u3-1) (1-t)^u4 and Q(s,t) = u1 s^u2 (1-s)^(u1-1) t^u3 (1-t)^u4 on the open square. Direct differentiation gives
 dP/ds = u2 u3 A - u1 u3 X, dQ/dt = u1 u3 X - u1 u4 C.
The cross terms +-u1u3 X cancel, so dP/ds + dQ/dt = u2 u3 A - u1 u4 C. Both derivatives lie in L^1(T) because A, X and C do.
(a) Fix t in (0,1). On s in [0,1-t] we have 1-s >= t > 0, so P(.,t) is continuous on [0,1-t] and C^1 on (0,1-t]. It vanishes at s=0 because u2 > 0, and dP/ds is integrable (it is at most a constant times s^(u2-1) near 0). Applying the fundamental theorem of calculus on [eps,1-t] and letting eps -> 0 gives integral_0^(1-t) dP/ds ds = P(1-t,t) = u3 t^(u1+u3-1) (1-t)^(u2+u4). By Fubini, integral over T of dP/ds = u3 Beta(u1+u3, 1+u2+u4).
(b) Fix s in (0,1). On t in [0,1-s] we have 1-t >= s > 0. The same argument, using Q(s,0) = 0 because u3 > 0, gives integral_0^(1-s) dQ/dt dt = Q(s,1-s) = u1 s^(u2+u4) (1-s)^(u1+u3-1). By Fubini, integral over T of dQ/dt = u1 Beta(1+u2+u4, u1+u3).
Both Beta integrals converge since u1+u3 > 0 and 1+u2+u4 > 0. No corner terms occur, because the argument uses only one-dimensional slices and never takes a limit along an arc.
Adding (a) and (b) gives (u1+u3) Beta(u1+u3, 1+u2+u4) = G(1+u1+u3)G(1+u2+u4)/G(1+u1+u2+u3+u4) = G0, using (u1+u3)G(u1+u3) = G(1+u1+u3). The left side is integral over T of (u2u3 A - u1u4 C) = F23 - F32. Hence F23 - F32 = G0.

Step 3 (the only use of the plane hypothesis): I_up = F32 on Pi1 or Pi2.
The map (s,t) -> (1-t,1-s) is a linear involution of the unit square with |det| = 1, and s'+t' = 2-(s+t). It therefore carries {s+t >= 1} onto T. Substituting s = 1-t' and t = 1-s' gives
 I_up = u2 u3 * integral over T of s^u4 (1-s)^(u3-1) t^u1 (1-t)^(u2-1).
On Pi2 (u1=u3, u2=u4): this is u1 u2 * integral over T of s^u2 (1-s)^(u1-1) t^u1 (1-t)^(u2-1). That equals F32 exactly, since u1u4 = u1u2, u3 = u1 and u4 = u2.
On Pi1 (u1=u2=a, u3=u4=b): this is ab * integral over T of s^b (1-s)^(b-1) t^a (1-t)^(a-1). F32 equals ab * integral over T of s^a (1-s)^(a-1) t^b (1-t)^(b-1). T is symmetric under s <-> t, so the two integrals agree.

Conclusion. On Pi1 or Pi2, F23 = V - I_up = V - F32 = V - (F23 - G0). Hence 2 F23 = V + G0, that is F23 = (1/2)(G0' + G0). QED.

Scope. Proved for real u_i > 0 on Pi1 or Pi2. Step 2 holds for every u in (0,inf)^4, and the plane hypothesis enters only through Step 3. The formula is false off the planes. For example F23(2,1,3,1) = 5/84 while the formula gives 11/168, and F23(1,2,3,4) = 2/315 while the formula gives 1/140. The proof does not show that Pi1 and Pi2 are the only loci where the formula holds.

## Analytic continuation

Established: real u_i > 0 on Pi1 or Pi2. The identical argument (all estimates depend only on |s^(u-1)| = s^(Re u - 1)) covers complex u with Re u_i > 0 on Pi1 or Pi2, with the plane relations read as complex equalities. The rigor lens confirmed this numerically to about 1e-23 at three complex plane points. Both sides are holomorphic there, all Gamma arguments have real part > 1, and there are no poles. NOT established: any continuation to Re u_i <= 0. The integral F23 diverges at the edges s=0 (Re u2 <= 0) and t=0 (Re u3 <= 0). The right-hand side is meromorphic on C^2, but the meromorphic continuation of F23 itself was not verified, and nothing beyond Re u_i > 0 should be claimed. Off the planes the formula fails, and nobody proved that Pi1 and Pi2 are the only planes where it holds.

## Residual doubts (judge)
- Exclusivity: nobody proved that Pi1 and Pi2 are the only loci where the formula holds. It fails numerically at about 12 off-plane points (blind lens, 36 digits) and in my exact check at (2,1,3,1), (1,2,3,4) and (3,1,2,2). Do not write 'if and only if'. A palindromic family u1=u4, u2=u3 also fails.
- The name F23 in this proof is the triangle integral in the u_i variables, tied to Vol(B)/16. The citations use F^(23) of Schlotterer-Stieberger, which is a unit-square integral in Mandelstam variables s_ij. No referee established that the two coincide under any dictionary u_i <-> s_ij. Do not claim identity, or that the 5-point cyclic relation (6.12) of Broedel-Schlotterer-Stieberger applies to this paper's F23 as a theorem, without a separate argument.
- F23 = Vol(B)/16 was verified by the rigor lens by hand and by Monte Carlo (about 1.4 sigma on the planes; off-plane MC matches F23 and not the claim formula). The blind lens did not check it. I re-derived the change of variables above, and it is elementary, but it is outside the claim.
- Complex u and analytic continuation: only Re u_i > 0 on the planes is covered (see analytic_continuation_domain). The meromorphic continuation of F23 beyond that domain was not verified.
- The numeric evidence is not a proof, and every statement above rests on the algebraic proof. Caution for anyone rerunning the numerics: mpmath hyp3f2 at z=1 is unreliable for large non-integer parameters. It was wrong at 12 of 16 points at dps 60, and not because of the claim. Direct term-by-term summation or quadrature at high precision agrees with the claim. mp.quad also has an effectively absolute tolerance, which makes relative errors look poor when F23 is below about 1e-40. Plain tanh-sinh on t^(u3-1) endpoint singularities produced spurious residuals until the singularities were removed by substitution.
- The blind lens gave a second, independent proof route (F_T - F_T' as a sgn(r-t) double integral, with separate integrations by parts on Pi1 and Pi2), checked numerically at 36 digits. It was not refereed by anyone else, and I did not include it. The proof above does not depend on it.
- Citations: all equation numbers are arXiv-version numbers (1205.1516 v3, 1304.7267 v2, 1908.10830 v2, 1304.7304 v2). Published-journal numbering was not checked for any of them. The coefficient 1/5 of Q_8 was read, not recomputed, and is basis-dependent. The Mafra-Schlotterer dropout of zeta_{3,5}, zeta_{3,7} and zeta_{3,5,3} was confirmed only as a stated remark, with no proof or checked order given. The first names Oliver and Stephan in the Schlotterer-Stieberger entry come from companion records, not the J. Phys. A record itself.
- Housekeeping: the referees left scratch files in <tmp> (including its sonnet_verify and scripts subdirectories). They are safe to delete once the parent is finished with them, and I did not touch them. My own check used only an inline python -c command and wrote no files. No repo file was edited.

Planner note on the identification doubt: with (s12,s23,s34,s45,s24) = (u3,0,u2,u1+u3,u4), the printed closed form of Schlotterer-Stieberger eq. (3.11) for F^(23) becomes exactly G*3F2(-u4,u3,1+u1+u3;1+u3,1+u1+u2+u3;1) = F23 (Gamma prefactor G(1+u2)G(1+u1+u3)/G(1+u1+u2+u3), same 3F2 parameters); the citation referee checked (3.11) numerically against its Euler integrals. So F23 = F^(23)|_{s23=0} via the printed formula.
