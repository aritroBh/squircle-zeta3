code = r'''
import mpmath as mp, pickle, time
mp.mp.dps = 620
t0 = time.time()
D1 = '../DOD1_dodecagon/'
r3 = mp.sqrt(3)
X = {1: 1-2*r3/3, 2: 1-r3, 5: 1-r3/3, 6: 1-r3/2}
IDX = [1,2,5,6]
atoms = pickle.load(open(D1+'atoms620.pkl','rb'))
T3a = atoms['T3a']; Lam620 = atoms['Lam']

# gamma_{i,j} = x_j/(x_j - x_i), delta_{i,j} = gamma_{i,j}(1-x_i)
gam = {}; dl = {}
for i in IDX:
    for j in IDX:
        if i==j: continue
        gam[(i,j)] = X[j]/(X[j]-X[i])
        dl[(i,j)]  = gam[(i,j)]*(1-X[i])

Li2R = lambda z: mp.re(mp.polylog(2, mp.mpc(z)))
Li3R = lambda z: mp.re(mp.polylog(3, mp.mpc(z)))
def Creal(a,b):
    return mp.quad(lambda t: mp.log(abs(1-a*t))*mp.log(abs(1-b*t))/t, [0,1])
def SwpR(a,b):
    # real swap part of C(a,b): = Re Li_{2,1}(b/a, a)
    return Creal(a,b) - mp.log(abs(1-a))*(Li2R(b) - Li2R(b/a)) - Li3R(b)
def Li21ser(x,y):
    # convergent series for |x|<1,|y|<1
    tot=mp.mpf(0); ps=mp.mpf(0); xp=mp.mpf(1); yp=mp.mpf(1)
    for m in range(1,200000):
        xp*=x
        if m>=2: tot += xp*ps/m**2
        yp*=y; ps += yp/m
        if m>50 and abs(xp*ps)<mp.mpf(10)**(-660): break
    return tot

triples = sorted(T3a.keys())
SwG = {}; SwD = {}
worst_ser = mp.mpf(0)
for (i,j,k) in triples:
    a1,b1 = gam[(i,j)], gam[(i,k)]
    a2,b2 = dl[(i,j)], dl[(i,k)]
    sg = SwpR(a1,b1); sd = SwpR(a2,b2)
    SwG[(i,j,k)] = sg; SwD[(i,j,k)] = sd
    # cross-check vs series when both args in (-1,1)
    for (a,b,sw,tag) in [(a1,b1,sg,'gam'),(a2,b2,sd,'del')]:
        if abs(a)<1 and abs(b/a)<1:
            d = abs(sw - Li21ser(b/a,a))
            worst_ser = max(worst_ser,d)
    print((i,j,k), "SwG=%s SwD=%s" % (mp.nstr(sg,12), mp.nstr(sd,12)), "t=%ds"%(time.time()-t0))
print("worst |SwpR - series| over convergent cases:", mp.nstr(worst_ser,4))
pickle.dump({'gam':gam,'dl':dl,'SwG':SwG,'SwD':SwD}, open('gdswaps620.pkl','wb'))
print("saved gdswaps620.pkl, elapsed", time.time()-t0)
'''
open('step_03.py','w').write(code)
exec(code)
