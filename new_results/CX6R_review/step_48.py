code = r'''import mpmath as mp
import cx6r_common as C
print("CX6R step_16d: sec.3 generating-function identity, well-conditioned reference")
mp.mp.dps = 80
def B(z): return 1+mp.exp(-z)-mp.sqrt(1+mp.exp(-2*z))
def Cc(z): return mp.exp(z)*B(z)
def M_form2(eps, level=6):
    def f(x):
        return mp.exp(-x)*(1+mp.exp(-x))**(-1-eps)*Cc((1-eps)*x)
    return -C.glint(f, [[0,1],[1,4],[4,12],[12,40],[40,150]], level)
def Phi(eps, s):
    return C.glint(lambda u: mp.exp(-(s*(1-eps)+1)*u)*(1+mp.exp(-u))**(-1-eps),
                   [[0,1],[1,3],[3,9],[9,30]], 6)
def J(eps): return (1-mp.mpf(2)**(-eps))/eps if eps != 0 else mp.log(2)
for eps in [mp.mpf('0.05'), mp.mpf('-0.1'), mp.mpf('0.2'), mp.mpf('-0.24'), mp.mpf('0.001')]:
    I = C.glint(lambda r: C.H(mp.mpf('0.5')+r*C.w4)*C.K(mp.mpf('0.5')+r*C.w4)*Phi(eps, mp.mpf('0.5')+r*C.w4)*C.w4, C.GAMMA_SEGS, 6)
    val = -J(eps) + 2/mp.pi*mp.im(I)
    print(f"  eps={mp.nstr(eps,6)}: |M(x-form) - (-J + 2/pi Im int_Gamma)| = {mp.nstr(abs(val-M_form2(eps)),5)}")
'''
open('step_16d.py','w').write(code)
exec(compile(code,'step_16d.py','exec'))
