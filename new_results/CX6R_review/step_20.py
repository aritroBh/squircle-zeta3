code = r'''import mpmath as mp, re
mp.mp.dps = 80
w = mp.expjpi(mp.mpf(1)/4)
def H(s): return mp.sqrt(mp.pi)/4*mp.gamma(s/2)/mp.gamma(s/2+mp.mpf(3)/2)
def K(s): return 1/(2*mp.cos(mp.pi*s/2))
def Pk(k, s):
    # (-1)^k/k! int_0^1 t^s/(1+t) (s log t + log(1+t))^k dt
    def f(t):
        if t == 0: return mp.mpf(0)
        return t**s/(1+t)*(s*mp.log(t)+mp.log1p(t))**k
    pts = [mp.mpf(0)]+[mp.mpf(10)**(-e) for e in (6,4,2)][::-1]+[mp.mpf('0.5'), mp.mpf(1)]
    return (-1)**k/mp.factorial(k)*mp.quad(f, pts)
def Jk(k): return (-1)**k*mp.log(2)**(k+1)/mp.factorial(k+1)
# table M_k
tab = {}
for line in open('<project>/research/08_theory/D13f_sweep/mk_deg9.tsv'):
    p = line.rstrip('\n').split('\t')
    if p[0].isdigit():
        tab[int(p[0])] = mp.mpf(re.match(r'\[([-0-9.eE+]+)', p[1]).group(1))
print("CX6R step_05: eq (7) M_k = -J_k + (2/pi) Im int_Gamma H K P_k ds, vs table")
for k in [0,1,2,3,4,5,6,7,8,20,50]:
    def g(r):
        s = mp.mpf('0.5')+r*w
        return H(s)*K(s)*Pk(k,s)*w
    I = mp.quad(g, [0, 1, 4, 16, 64, mp.inf])
    val = -Jk(k) + 2/mp.pi*mp.im(I)
    err = abs(val-tab[k])
    rel = err/max(abs(tab[k]),mp.mpf('1e-60'))
    print(f"k={k:3d}  table={mp.nstr(tab[k],25)}")
    print(f"       eq(7)={mp.nstr(val,25)}   abs err={mp.nstr(err,5)}  rel={mp.nstr(rel,5)}")
'''
open('step_05.py','w').write(code)
exec(compile(code,'step_05.py','exec'))
