
import os
os.chdir('<project>/research/08_theory/COR8_depth_two_any_angle')
res = open('result.md','w')
res.write('''# COR8: state of play (interrupted run)

## Established (PROVED at symbol level, exact rational arithmetic)
1. delta_22(Omega) = (1-c^-1)([c]^[1-c])^([1-c]^[1+c]) + (-1-c^-1)([c]^[1+c])^([1-c]^[1+c]) (recomputed, = COR4).
2. TRUE constant-coefficient symbols of Li_{3,1}(x,y) (18 terms) and Li_{2,2}(x,y) (15 terms) derived and validated
   (y=+-1 reductions to H_{0,0,1,+-1}, H_{0,1,0,+-1}; stuffle; Li4 and shuffles killed by d22). Engine: cor8_sym.py, step_01.py (all PASS).
3. delta_22 MATCH (exact): d22(Omega) = d22( ((1-c)/c) Li_{2,2}(1/2,1-c) - ((1+c)/c) Li_{2,2}(1/2,1+c) ).
   Two depth-two terms, both with x=1/2, y=1∓c; series convergent on (0,1). (step_02)
4. BUG FOUND in COR4 sunit_args.json: 20 of 40 [1-r] vectors have the [1-c] exponent sign flipped
   (step13 factor_vec treated (c-1)^e as -e[1-c]; [c-1]=[1-c] so it must be +e). Verified by independent factorisation.
5. With corrected vectors: exact symbol identity S(Omega - DT) = sum a_i(c) S(Li4(g_i)) + sum p_k(c) S(products),
   residual EXACTLY 0 in all 243 components (81 words x 3 exponents). ~24 Li4 + ~30 product terms, small rational coeffs.
   (solve81_v2.pkl, full_sols_v2.pkl)

## NOT established / unfinished
- The function-level numeric verification of the 64-term identity FAILS (residuals ~30 at c=0.5, ~194 at c=0.13).
  Symbol residual is exactly 0, Omega evaluation independently verified (matches M3 decomposition to 40 digits at c=1/2),
  so the failure is either a remaining bug in the numeric RHS evaluator or a deeper convention issue.
- Steps 3-5 of the brief (constants, corners, final PROVED/NUMERICAL identity) not completed.
''')
res.close()
open('progress.log','a').write('step_59: result.md written; run interrupted with function-level check failing\n')
print(open('result.md').read())
