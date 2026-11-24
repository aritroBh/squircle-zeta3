import os
os.chdir('<project>/research/08_theory/RST1_rest_exponent')
p = open('proof.md').read()
p = p.replace(r"""is $h(k)\ge0$ with
$h(k)=\sqrt{\pi(k/2)(k/2+1)}-\sqrt{14\pi(k-13)}-\frac{a(k/2-12)}{2\sqrt2}$:
$h(99)=12.390>0$ (arb) and $h$ is increasing on $[99,300)$ (grid) and for $k\ge300$ by
$h(k)\ge\frac{\sqrt\pi}{2}k-\sqrt{14\pi k}-\frac{ak}{4\sqrt2}\ge k\bigl(\frac{\sqrt\pi}{2}-\frac{a}{4\sqrt2}\bigr)-\sqrt{14\pi k}>0$
for $k\ge300$.""",
r"""is $\ell(13)\ge\ell(k/2)$ with $\ell(i)=\frac{a(i+1)}{2\sqrt2}-\sqrt{\pi(k-i)(i+1)}$,
i.e. $h(k)\ge0$ with
$h(k)=\sqrt{\pi(k/2)(k/2+1)}-\sqrt{14\pi(k-13)}-\frac{a(k/2-13)}{2\sqrt2}$:
$h(99)=12.78>0$ (arb); $h(k)>0$ for $k\in[99,400]$ (direct arb check of
$Q_{13}\ge Q_{k/2}$) and for $k>400$ by
$h(k)\ge\frac{\sqrt\pi}{2}k-\sqrt{14\pi k}-\frac{ak}{4\sqrt2}
\ge0.69\,k-\sqrt{14\pi k}>0$ for $k\ge300$.""")
# also clarify the e^{a/2} in the mid2 lump (13)
p = p.replace(r"""\[
\sum_{i=13}^{k/2}\frac{2^i}{i!}F^{\rm box}_{H,a}(k-i,i+1)
\le\frac{2\sqrt2}{14}\Bigl(2\sqrt{\tfrac{k(k+1)}{\pi}}\,Q_{13}
+\frac1a e^{-\sqrt{14\pi(k-13)}}\Bigr)\sum_{i\ge13}\frac{2^i}{i!}.\tag{13}
\]""",
r"""(the box prefactor is $e^{ac_0(i)/2}=e^{a/2}e^{\frac{a(i+1)}{2\sqrt2}}$, so the
prefactor–exponent product is $e^{a/2}Q_i\le e^{a/2}Q_{13}$ by R5, the factor $e^{a/2}$
being $i$-independent; the displayed $Q_{13}$ below already includes it):
\[
\sum_{i=13}^{k/2}\frac{2^i}{i!}F^{\rm box}_{H,a}(k-i,i+1)
\le\frac{2\sqrt2}{14}\Bigl(2\sqrt{\tfrac{k(k+1)}{\pi}}\,e^{\frac{a(1+14/\sqrt2)}{2}}
e^{-\sqrt{14\pi(k-13)}}+\frac1a e^{-\sqrt{14\pi(k-13)}}\Bigr)\sum_{i\ge13}\frac{2^i}{i!}.\tag{13}
\]""")
open('proof.md','w').write(p)
print("proof.md updated")
# fix verify_rst1.py hfun to -13
v = open('verify_rst1.py').read()
v = v.replace("return (PI*(k/2)*(k/2+1)).sqrt()-(14*PI*(k-13)).sqrt()-a*(k/2-12)/(2*SQ2)",
              "return (PI*(k/2)*(k/2+1)).sqrt()-(14*PI*(k-13)).sqrt()-a*(k/2-13)/(2*SQ2)")
v = v.replace('log(f"CERTIFIED: Q_i<=Q_2 on [2,k/2] ({bad} violations); Q_13>=Q_(k/2) on [99,400] ({bad2}); h(99)={hfun(99).str(8)}>0")',
              'log(f"CERTIFIED: Q_i<=Q_2 on [2,k/2] ({bad} violations); Q_13>=Q_(k/2) on [99,400] ({bad2}); h(99)={hfun(99).str(8)}>0, h(k)>0 for k>=99")')
open('verify_rst1.py','w').write(v)
print("verify_rst1.py updated")
import subprocess
r = subprocess.run(['bash','-c','''cd <project>/research/08_theory/RST1_rest_exponent && <project>/.venv/bin/python -c "
from flint import arb, ctx
ctx.dps=50
PI=arb.pi(); SQ2=arb(2).sqrt(); a=PI/(2*SQ2)
def h(k):
    k=arb(k)
    return (PI*(k/2)*(k/2+1)).sqrt()-(14*PI*(k-13)).sqrt()-a*(k/2-13)/(2*SQ2)
print('h(99)=',h(99).str(8), ' h(400)=',h(400).str(8))
"'''], capture_output=True, text=True, timeout=300)
print(r.stdout, r.stderr[-300:])
