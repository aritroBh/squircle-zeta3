
import os
os.chdir('<project>/research/08_theory/COR8_depth_two_any_angle')
res = open('result.md','w')
res.write('''# COR8: state of play at end of run

## Established (PROVED at symbol level; exact rational arithmetic throughout)
1. **delta_22(Omega)** recomputed from the 74-term HPL fit = COR4's value:
   (1-c^-1)([c]^[1-c])^([1-c]^[1+c]) + (-1-c^-1)([c]^[1+c])^([1-c]^[1+c]). (step_01)
2. **TRUE constant-coefficient symbols** of Li_{3,1}(x,y) (18 terms) and Li_{2,2}(x,y) (15 terms), derived from the
   differential recursion with 1-forms converted to exact dlog form. Validated: y=1 -> H_{0,0,1,1}/H_{0,1,0,1};
   y=-1 -> -H_{0,0,1,-1}/-H_{0,1,0,-1}; stuffle Li11(x,y)+Li11(y,x)=Li1x Li1y - Li2(xy); d22(Li4)=d22(shuffles)=0;
   MPL series = genuine double sums (0.0 residual at 60 dps). (cor8_sym.py, step_01/61)
3. **delta_22 answered with TWO depth-two terms** (exact, two equivalent forms):
   DT_A = ((1-c)/c) Li_{2,2}(1/2,1-c) - ((1+c)/c) Li_{2,2}(1/2,1+c)   [used in the solve]
   DT_B = ((1-c)/c) Li_{2,2}(1/2,1-c) + ((1+c)/c) Li_{3,1}(1/2,1+c)
   Both give d22(Omega - DT) = 0 exactly (3-letter symbols). The DGR 3-term variant also analysed. (step_02, step_60)
4. **COR4 sunit_args.json BUG**: 20 of 40 [1-r] vectors have the [1-c] exponent sign flipped (step13 factor_vec
   treated (c-1)^e as -e[1-c] instead of +e). Verified by independent factorisation; corrected tables rebuilt.
5. With [2] dropped, an exact symbol identity S(Omega - DT) = sum a_i(c) S(Li4(g_i)) + sum p_k(c) S(products)
   exists, residual EXACTLY 0 in all 243 components (~19 Li4 + ~33 products, small rational coefficients).
   Its function-level numeric check FAILS (residual ~30 at c=1/2) because symbol-invisible terms are not yet
   included (Step 4: log2*(w3), zeta(2)*(w2), zeta(3)*(w1), zeta(4)). (full_sols_v3.pkl)
6. Numeric infrastructure validated: Omega(1/2) matches M3 decomposition to 40 digits; Re(polylog)=Re(inversion)
   at x=2,4 to 60 digits; Puiseux-series engine cor8_series.py validated for all 24 args to 40+ digits.
   Residual series D computed: c^-1(zeta(2) log^2 c + const) + regular log^0/log^1 terms.
7. Level-2 constant reductions verified at 120 digits:
   Li_{2,2}(1/2,1) = (3 zeta(4) + 12 zeta(3) log 2 - 12 zeta(2) log^2 2 + 2 log^4 2)/48
   Li_{3,1}(1/2,1) = zeta(4)/8 - zeta(3) log 2 / 8 + Li_4(1/2)/24

## In progress / NOT established
- With [2] kept (to make log2-terms visible), DT leaves a [2]-residual delta_22 with 4 components. The full
  [2]-kept ansatz (38 Li4 + 1212 products + 209 log2-terms + 448 depth-two symbols, rank 121 < 256) does NOT yet
  contain the residual: 2 violated left-null constraints at E=-1 (patterns m2cm/mcm2/mcmp/mcpm and
  2mcm/2pcm/m2cm/pcm2/pcmp/pcpm) whose parent functions were NOT found among all 224 S-unit depth-two pairs
  (enumerated and pairing-checked). A pair ((1-c)/2, 1/(1+c)) that seemed to fit had 1-xy = (1+3c)/(2(1+c)),
  outside the alphabet (spurious).
- Hence: the full function identity, the corner checks (c=1/sqrt2 vs thm:W4R, c=1/2 vs COR5 1095 digits, c=1/3),
  simplification of the Li4 tail, verify_cor8.py and proof.md are UNFINISHED.
- No function-level identity is claimed.
''')
res.close()
open('progress.log','a').write('final: run ended; result.md updated. [2]-kept solve lacks 2 constraints; unfinished.\n')
print("result.md written")
