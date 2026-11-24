
fix4 = r'''
import mpmath as mp
mp.mp.dps = 80
def Vnum(etaval):
    q = 1+mp.mpf(etaval)
    f = lambda v: (2/q)*(v*(1-v))**(1/q-1)*(v**(q-1)+(1-v)**(q-1))**(1/q)
    return 2*mp.quad(f, [mp.mpf('1e-14'), mp.mpf('1e-8'), mp.mpf('1e-4'), mp.mpf('0.01'), mp.mpf('0.5')])
print("Duality with careful (symmetry-doubled) quadrature:")
for ev in ['0.25','0.1','0.05','1','3']:
    a = Vnum(ev); b = Vnum(mp.mpf(1)/mp.mpf(ev))
    print(f"eta={ev}: V={mp.nstr(a,25)}  V(1/eta)={mp.nstr(b,25)}  diff={mp.nstr(abs(a-b),3)}")
'''
open('fix4.py','w').write(fix4)
import subprocess
r = subprocess.run(['<project>/.venv/bin/python','fix4.py'],capture_output=True,text=True,timeout=900)
print(r.stdout); print(r.stderr[-300:])
