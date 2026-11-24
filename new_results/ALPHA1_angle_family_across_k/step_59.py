open('progress.log','a').write('step_20: T(alpha) = pi^2/12 - Li2(c) - Li2(-c) = pi^2/12 - (1/2)Li2(c^2), 16 angles 88 digits + 3 holdouts; reduces to squircle T\n')
# quick sanity of the Li2 duplication and squircle specialization
import subprocess
PY = '<project>/.venv/bin/python'
chk = r'''
import mpmath as mp
mp.mp.dps = 50
c = 1/mp.sqrt(2)
lhs = mp.polylog(2, c) + mp.polylog(2, -c)
rhs = mp.polylog(2, c*c)/2
print("dup:", mp.nstr(abs(lhs-rhs), 3))
T = mp.pi**2/12 - lhs
print("T(pi/2) =", mp.nstr(T, 30), " target:", mp.nstr((mp.pi**2 + 6*mp.log(2)**2)/24, 30))
'''
r = subprocess.run(['nice','-n','10',PY,'-c',chk], capture_output=True, text=True)
print(r.stdout)

src = open('result.md').read()
src += r'''

## 4. Weight two (Task B, part 1) — DONE

THEOREM (EXACT for the kernel at 4 rational angles, k<=14; NUMERICAL 76-88 digits for the value):
For every k >= 2 the degree-1 kernel is the UNIVERSAL (angle-independent) rational function
  K^{(k)}_{10}(w) = 1/(2w) + 1/(2(w-1)) - 1/(w+1) - 1/(2(w+r)),
i.e. the squircle's R_w2 of eq:E mapped by v = a w.  Hence W_2(k;alpha) = T(alpha) is constant in k for k>=2
(checked EXACTLY at c = 1/2, 3/5, 1/3, 4/5 for every k = 2..14: residues identical at all k, all angles).
Equivalently the epsilon-series identity (1-eps)(W_10 - eps R_10) = eps^2 R_w2 holds at every tested angle.

VALUE (NUMERICAL, 16 angles to 88 digits + 3 holdouts to 76 digits, step_20; verified also vs corner_series):
  T(alpha) = pi^2/12 - Li_2(c) - Li_2(-c) = pi^2/12 - (1/2) Li_2(c^2),   c = cos(alpha/2).
CONTROL: at c = 1/sqrt2, Li2(1/sqrt2)+Li2(-1/sqrt2) = (1/2)Li2(1/2) gives T = (pi^2 + 6 log^2 2)/24, the squircle T.
'''
open('result.md','w').write(src)
print("recorded")
