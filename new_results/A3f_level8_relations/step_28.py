src=open('step_08_weight2_certificate.py').read()
src=src.replace("v.get(cols[q], F(0)) * F(rr[i, j])", "v.get(cols[q], F(0)) * F(str(rr[i, j]))")
open('step_08_weight2_certificate.py','w').write(src)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python','-B','step_08_weight2_certificate.py'],
                   capture_output=True, text=True, timeout=2400)
print(r.stdout[-2500:]); print('STDERR:', r.stderr[-1500:])
open('progress.log','a').write(
 'step_05/06: weight-2 values at 325 digits; LLL found missing relation [0,-6,-1,7,3,-9,-3,21,-5,5] (residual 1.5e-329).\n'
 'step_07: 20 Rogers five-term + 5 reflection instances on real cross-ratios; rank 54->55; quotient 9. All rows numeric zero (worst 1.2e-58).\n'
)
