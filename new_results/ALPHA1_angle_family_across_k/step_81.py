open('progress.log','a').write('step_30: W3 slope = s(alpha)*J(c), s=(1-c^2)/(2c^2); J exact HPL (6 terms, verified sq squircle); classical closed form verified 61 digits x6 angles\n')
src = open('result.md').read()
src += r'''

## 5. Weight three (Task B, part 2) — DONE

THEOREM (EXACT for the kernel law at 4 rational angles, k = 3..14; NUMERICAL 61 digits for the value):
For k >= 3 the degree-2 kernels are affine in k, with slope kernels UNIVERSAL times s(alpha):
  S_20 = dK^{(k)}_{20}/dk = s(alpha) (1/(w+1) - 1/(w+r)),   S_11 = dK^{(k)}_{11}/dk = s(alpha) (1/w - 2/(w+1)),
  s(alpha) = (1-c^2)/(2c^2) = (1/2) tan^2(alpha/2)   [EXACT at c=1/2,3/5,1/3,4/5 for every k=3..14;
  at c=1/sqrt2, s=1/2, these are the squircle slope kernels S_20, S_11 of thm:w3k].
So W_3(k;alpha) = W_3(3;alpha) + (k-3) slope(alpha), with slope(alpha) = s(alpha) J(c), and (NUMERICAL, 61 digits
at 6 angles incl. squircle control; exact HPL form step_30, verified vs quadrature):
  J(c) = log(1-c^2) [Li_2(c^2) - zeta(2)] + 2 Li_3(1-c^2) - 2 zeta(3) + log(c^2) log^2(1-c^2)
       = 2 T(alpha) log(1-c^2) - I(c^2),   I(x) = int_0^x log^2(1-u)/u du.
Equivalently J = zeta(2) log(1-c^2) + int_0^{c^2} Li_2(u)/(1-u) du as an exact 6-term HPL in c
(zeta(2)(H_{-1}-H_1) + 2(H_{-1,0,-1}-H_{-1,0,1}-H_{1,0,-1}+H_{1,0,1})).
CONTROL: at c=1/sqrt2, slope = -(zeta(3)/8 + pi^2 log2/24 + log^3 2/12), the squircle value (61 digits).
So the slope IS a closed form in c (Li_3(1-c^2), Li_2(c^2), logs of c^2, 1-c^2).
'''
open('result.md','w').write(src)
print("recorded weight three")
