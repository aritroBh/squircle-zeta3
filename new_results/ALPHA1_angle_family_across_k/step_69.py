step28 = r'''"""step_28.py -- J(r) identification: high-precision PSLQ at c=1/2 (r=1/3), verify at c=3/5, 1/3."""
import sys
sys.path.insert(0, '.')
import mpmath as mp
mp.mp.dps = 150

def Jval(c):
    r = (1-c)/(1+c)
    def lt(w): return mp.log((w-r)*(w+1)/(2*(1-r)*w))
    def lmt(w): return mp.log((1-w)*(w+r)/(2*(1-r)*w))
    def f(w):
        return (1/(w+1) - 1/(w+r))*lt(w)**2 + (1/w - 2/(w+1))*lt(w)*lmt(w)
    return mp.quad(f, [r, mp.sqrt(r), 1])

J = Jval(mp.mpf(1)/2)
print("J(1/3) =", mp.nstr(J, 60))
L2, L3 = mp.log(2), mp.log(3)
A3 = lambda x: mp.polylog(3, x)
B = [mp.zeta(3), A3(mp.mpf(1)/3), A3(-mp.mpf(1)/3), A3(mp.mpf(2)/3), A3(mp.mpf(1)/2), A3(-mp.mpf(1)/2),
     A3(mp.mpf(1)/4), A3(mp.mpf(3)/4), A3(-mp.mpf(2)/3), A3(-mp.mpf(1)/4), A3(mp.mpf(1)/9),
     mp.pi**2*L2, mp.pi**2*L3, L2**3, L3**3, L2**2*L3, L2*L3**2,
     mp.polylog(2, mp.mpf(1)/3)*L3, mp.polylog(2, mp.mpf(1)/3)*L2, mp.polylog(2, -mp.mpf(1)/3)*L2]
names = ['zeta3', 'Li3(1/3)', 'Li3(-1/3)', 'Li3(2/3)', 'Li3(1/2)', 'Li3(-1/2)', 'Li3(1/4)', 'Li3(3/4)',
         'Li3(-2/3)', 'Li3(-1/4)', 'Li3(1/9)', 'pi^2 l2', 'pi^2 l3', 'l2^3', 'l3^3', 'l2^2 l3', 'l2 l3^2',
         'Li2(1/3) l3', 'Li2(1/3) l2', 'Li2(-1/3) l2']
v = [J] + B
rel = mp.pslq(v, tol=mp.mpf(10)**(-135), maxcoeff=10**8, maxsteps=5000)
if rel:
    terms = [(cf, names[i-1]) for i, cf in enumerate(rel) if cf and i >= 1]
    print("relation (rel[0] =", rel[0], "):", terms)
    if rel[0] != 0:
        scale = -rel[0]
        def Jpred(c):
            r = (1-c)/(1+c)
            bmap = {'zeta3': mp.zeta(3)}
            pred = mp.mpf(0)
            for cf, nm in terms:
                pass
            return pred
        # verify at c=3/5, 1/3 by direct evaluation of the found combination
        for c in (mp.mpf(3)/5, mp.mpf(1)/3):
            r = (1-c)/(1+c)
            lc, l1c, lmc, l2v = mp.log(c), mp.log(1+c), mp.log(1-c), mp.log(2)
            lr, l1r, lmr = mp.log(r), mp.log(1+r), mp.log(1-r)
            A3v = lambda x: mp.polylog(3, x); A2v = lambda x: mp.polylog(2, x)
            vals = {'zeta3': mp.zeta(3), 'Li3(1/3)': A3v(r), 'Li3(-1/3)': A3v(-r), 'Li3(2/3)': A3v(1-r),
                    'Li3(1/2)': A3v(mp.mpf(1)/2), 'Li3(-1/2)': A3v(-mp.mpf(1)/2),
                    'Li3(1/4)': A3v((1-r)/2), 'Li3(3/4)': A3v((1+r)/2), 'Li3(-2/3)': A3v(r-1), 'Li3(-1/4)': A3v((r-1)/2),
                    'Li3(1/9)': A3v(r*r), 'pi^2 l2': mp.pi**2*l2v, 'pi^2 l3': mp.pi**2*mp.log(3),
                    'l2^3': l2v**3, 'l3^3': mp.log(3)**3, 'l2^2 l3': l2v**2*mp.log(3), 'l2 l3^2': l2v*mp.log(3)**2,
                    'Li2(1/3) l3': A2v(r)*lr, 'Li2(1/3) l2': A2v(r)*l2v, 'Li2(-1/3) l2': A2v(-r)*l2v}
            # NOTE: log names mapped to the c=1/2 logs; for verification reinterpret ln as the SAME r-based args
            pred = sum((mp.mpf(cf)/scale)*vals[nm] for cf, nm in terms)
            Jt = Jval(c)
            print(f"   verify c={mp.nstr(c,8)}: J={mp.nstr(Jt,30)} pred={mp.nstr(pred,30)} diff={mp.nstr(abs(Jt-pred),3)}")
else:
    print("no relation")
print("STEP28 DONE")
'''
open('step_28.py','w').write(step28)
import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run(['nice','-n','10',PY,'step_28.py'], capture_output=True, text=True, timeout=1700)
print(r.stdout)
print(r.stderr[-1500:] if r.returncode else "RC=0")
