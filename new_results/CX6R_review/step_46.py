code = r'''import mpmath as mp
import cx6r_common as C
print("CX6R step_16b: sec.3 generating-function identity, u-form Phi")
mp.mp.dps = 80
def M_orig(eps):
    def f(t):
        a = t**(eps-1); b = (1-t)**(eps-1)
        return a + b - a*mp.sqrt(1+(b/a)**2)
    return -C.glint(f, [[mp.mpf('1e-30'), mp.mpf('0.25')],[mp.mpf('0.25'), mp.mpf('0.5')]], 6)
def Phi(eps, s):
    # int_0^1 t^{s(1-eps)}(1+t)^{-1-eps} dt = int_0^inf e^{-(s(1-eps)+1)u}(1+e^{-u})^{-1-eps} du
    return C.glint(lambda u: mp.exp(-(s*(1-eps)+1)*u)*(1+mp.exp(-u))**(-1-eps),
                   [[0,1],[1,3],[3,9],[9,30]], 6)
def J(eps): return (1-mp.mpf(2)**(-eps))/eps if eps != 0 else mp.log(2)
for eps in [mp.mpf('0.05'), mp.mpf('-0.1'), mp.mpf('0.2'), mp.mpf('-0.24')]:
    I = C.glint(lambda r: C.H(mp.mpf('0.5')+r*C.w4)*C.K(mp.mpf('0.5')+r*C.w4)*Phi(eps, mp.mpf('0.5')+r*C.w4)*C.w4, C.GAMMA_SEGS, 6)
    val = -J(eps) + 2/mp.pi*mp.im(I)
    print(f"  eps={mp.nstr(eps,6)}: |M_orig - (-J + 2/pi Im int)| = {mp.nstr(abs(val-M_orig(eps)),5)}")
'''
open('step_16b.py','w').write(code)
exec(compile(code,'step_16b.py','exec'))
