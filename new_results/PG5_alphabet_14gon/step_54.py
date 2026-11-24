gold = r'''
"""PG5 step_06: independent gold P'''(1) of the regular 14-gon via polar perimeter.
K_p = {sum_i |<n_i,x>|^p <= 1}; r(th,p) = f(th,p)^(-1/p); P(p) = 2m int_chamber sqrt(r^2+r'^2).
9-point central finite difference for P'''(1). Fully independent of the E3n pipeline.
"""
import mpmath as mp, pickle
mp.mp.dps = 110
m = 7
def fsum(th, p):
    return mp.nsum(lambda i: abs(mp.cos(th - i*mp.pi/m))**p, [0, m-1])
def r_and_rp(th, p):
    f = fsum(th,p)
    fp = mp.nsum(lambda i: -mp.sin(th - i*mp.pi/m)*p*abs(mp.cos(th-i*mp.pi/m))**(p-1)*mp.sign(mp.cos(th-i*mp.pi/m)), [0,m-1])
    r = f**(-1/p)
    rp = -(1/p)*f**(-1/p-1)*fp
    return r, rp
def perim(p):
    # chamber (pi/2, pi/2+pi/m) times 2m ; inside chamber no kinks
    g = lambda th: (lambda rr: mp.sqrt(rr[0]**2+rr[1]**2))(r_and_rp(th,p))
    return 2*m*mp.quad(g, [mp.pi/2+mp.mpf('1e-6'), mp.pi/2+mp.pi/m-mp.mpf('1e-6')])
# central differences
h = mp.mpf('0.004')
vals = {}
for k in range(-4,5):
    vals[k] = perim(1+k*h)
    print("P(1%+dh) = %s"%(k, mp.nstr(vals[k],25)), flush=True)
# 9-point stencil for third derivative
c = { -4: mp.mpf(7)/240, -3: -mp.mpf(3)/10, -2: mp.mpf(169)/120, -1: -mp.mpf(61)/30, 0: mp.mpf(0),
      1: mp.mpf(61)/30, 2: -mp.mpf(169)/120, 3: mp.mpf(3)/10, 4: -mp.mpf(7)/240 }
# note: that's the second-derivative stencil; third derivative 8-point antisymmetric:
c3 = { -4: -mp.mpf(7)/240, -3: mp.mpf(3)/10, -2: -mp.mpf(169)/120, -1: mp.mpf(61)/30,
       1: -mp.mpf(61)/30, 2: mp.mpf(169)/120, 3: -mp.mpf(3)/10, 4: mp.mpf(7)/240 }
P3 = sum(c3[k]*vals[k] for k in c3)/h**3
print("P'''(1) gold =", mp.nstr(P3, 40))
# also 4-point and 6-point for error estimate
c3b = { -2: -mp.mpf(1)/2, -1: mp.mpf(1), 1: -mp.mpf(1), 2: mp.mpf(1)/2 }
P3b = sum(c3b[k]*vals[k] for k in c3b)/h**3
c3c = { -3: mp.mpf(1)/8, -2: -mp.mpf(1), -1: -mp.mpf(13)/8, 1: mp.mpf(13)/8, 2: mp.mpf(1), 3: -mp.mpf(1)/8 }
P3c = sum(c3c[k]*vals[k] for k in c3c)/h**3
print("P3 (4pt) =", mp.nstr(P3b,40))
print("P3 (6pt) =", mp.nstr(P3c,40))
print("spread:", mp.nstr(abs(P3-P3b),5), mp.nstr(abs(P3-P3c),5))
pickle.dump({'P3_gold':P3,'P3b':P3b,'P3c':P3c,'h':str(h),'vals':{k:str(v) for k,v in vals.items()}}, open('gold14.pkl','wb'))
open('progress.log','a').write('PG5 step_06: gold P3(1) 14-gon = %s (spread %s)\n'%(mp.nstr(P3,25),mp.nstr(abs(P3-P3c),3)))
'''
open('step_06_gold14.py','w').write(gold)
import subprocess, time
t0=time.time()
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_06_gold14.py'],capture_output=True,text=True,timeout=1200)
print("elapsed", time.time()-t0)
print(r.stdout); print(r.stderr[-500:] if r.stderr else '')
