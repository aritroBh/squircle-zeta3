import sys, os, time, pickle
from fractions import Fraction as Fr
work = '<project>/research/08_theory/TOW3_tower_depth6'
sys.path.insert(0, work); os.chdir(work)
os.nice(10)
import sympy as sp
import tow2lib as T
H = pickle.load(open('<project>/research/08_theory/TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb'))
def u(i,j,k,delta):
    return tuple(H[(i,j)][delta][k]) if (i,j) in H and delta in H[(i,j)] and k in H[(i,j)][delta] else (Fr(0),)*4
def gmul(x,y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def gdiv(x,y):
    den = y[0]**2+y[1]**2
    return ((x[0]*y[0]+x[1]*y[1])/den, (x[1]*y[0]-x[0]*y[1])/den)
def qsi_div(x, y):
    A=(x[0],x[2]); B=(x[1],x[3]); C=(y[0],y[2]); D=(y[1],y[3])
    C2m2D2 = (gmul(C,C)[0]-2*gmul(D,D)[0], gmul(C,C)[1]-2*gmul(D,D)[1])
    numA = (gmul(A,C)[0]-2*gmul(B,D)[0], gmul(A,C)[1]-2*gmul(B,D)[1])
    numB = (gmul(B,C)[0]-gmul(A,D)[0], gmul(B,C)[1]-gmul(A,D)[1])
    return (gdiv(numA, C2m2D2), gdiv(numB, C2m2D2))
a5 = sp.sympify(pickle.load(open('arho56_simplified.pkl','rb'))[5], locals={'d':T.d,'m':T.m,'I':sp.I})
a6 = sp.sympify(pickle.load(open('arho56_simplified.pkl','rb'))[6], locals={'d':T.d,'m':T.m,'I':sp.I})
dsym, msym = T.d, T.m
t0 = time.time()
bad = {5: [], 6: []}; tested = {5: 0, 6: 0}
for key in H:
    if key == '_meta': continue
    p, q = key; dd = p+q
    if 2*dd-4 < 6: continue
    e00 = u(p,q,dd,0)
    for delta, aa in [(5,a5),(6,a6)]:
        rv = sp.expand(aa.subs({dsym: dd, msym: p}))
        expected = (sp.Rational(sp.re(rv)), sp.Rational(sp.im(rv)))
        ratio = qsi_div(u(p,q,dd,delta), e00)
        tested[delta] += 1
        if ratio[1] != (Fr(0),Fr(0)) or ratio[0] != expected:
            bad[delta].append((key, ratio, expected))
print('rho(5): tested', tested[5], 'bad', len(bad[5]), '   rho(6): tested', tested[6], 'bad', len(bad[6]))
# also verify a_3, a_4 against same data as a sanity check of the method (TOW2 proved those)
from collections import Counter
def parts(rem, minpart, cur):
    if rem == 0:
        yield tuple(sorted(cur)); return
    for p in range(minpart, rem+1):
        yield from parts(rem-p, p, cur+[p])
L = pickle.load(open('L_parsed.pkl','rb'))
def bell(delta, L):
    tot = sp.Integer(0)
    for part in parts(delta, 1, []):
        C = Counter(part)
        coef = sp.Integer(1)
        for sz, mult in C.items():
            coef *= (L[sz]**mult) / sp.factorial(mult)
        tot += coef
    return tot
a4 = sp.cancel(bell(4, L))
bad4 = 0; t4 = 0
for key in H:
    if key == '_meta': continue
    p, q = key; dd = p+q
    if 2*dd-4 < 4: continue
    rv = sp.expand(a4.subs({dsym: dd, msym: p}))
    expected = (sp.Rational(sp.re(rv)), sp.Rational(sp.im(rv)))
    ratio = qsi_div(u(p,q,dd,4), u(p,q,dd,0))
    t4 += 1
    if ratio[1] != (Fr(0),Fr(0)) or ratio[0] != expected: bad4 += 1
print('sanity rho(4) (TOW2-proved): tested', t4, 'bad', bad4)
print('time', round(time.time()-t0,1))
with open('progress.log','a') as f:
    f.write('step_28-32: assembled a5,a6 via Bell from L_simplified (verified vs proved L(r)); rho(5),rho(6) vs data pending\n')
