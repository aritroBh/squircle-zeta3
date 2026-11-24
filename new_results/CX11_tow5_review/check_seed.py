"""Distinguish additive-constant seed blindness from arbitrary function perturbations."""
import sympy as s
d=s.symbols('d');Pi=s.prod(2*d-2*j-1 for j in range(1,6))
def step(h): return s.cancel((2*d-1)/(2*d-11)*h.subs(d,d+1)-h)
for seed,top,closing in [(s.Integer(1),None,0),(d**5/Pi,120/Pi,0),(d**6/Pi,None,720/Pi)]:
 h=seed
 for _ in range(5): h=step(h)
 if top is not None: assert s.cancel(h-top)==0
 h=step(h);assert s.cancel(h-closing)==0
 print('seed perturbation',seed,'sixth/closing response:',s.factor(h),flush=True)
print('PASS: constant shifts invisible; degree-five perturbation changes top; degree-six perturbation changes closing.',flush=True)
