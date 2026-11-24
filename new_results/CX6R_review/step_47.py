code = r'''import mpmath as mp
import cx6r_common as C
print("CX6R step_16c: which side limits? M_orig GL vs x-form reference")
mp.mp.dps = 80
def B(z): return 1+mp.exp(-z)-mp.sqrt(1+mp.exp(-2*z))
def Cc(z): return mp.exp(z)*B(z)
def M_form2(eps, level=6):  # well-conditioned x-form
    def f(x):
        return mp.exp(-x)*(1+mp.exp(-x))**(-1-eps)*Cc((1-eps)*x)
    return -C.glint(f, [[0,1],[1,4],[4,12],[12,40],[40,150]], level)
def M_orig(eps, level):
    def f(t):
        a = t**(eps-1); b = (1-t)**(eps-1)
        return a + b - a*mp.sqrt(1+(b/a)**2)
    return -C.glint(f, [[mp.mpf('1e-30'), mp.mpf('0.25')],[mp.mpf('0.25'), mp.mpf('0.5')]], level)
for eps in [mp.mpf('0.05'), mp.mpf('-0.1'), mp.mpf('0.2')]:
    ref = M_form2(eps)
    for lev in (5,6,7):
        d = abs(M_orig(eps, lev)-ref)
        print(f"eps={mp.nstr(eps,5)} lev={lev}: |M_orig_GL - x-form| = {mp.nstr(d,5)}")
'''
open('step_16c.py','w').write(code)
exec(compile(code,'step_16c.py','exec'))
