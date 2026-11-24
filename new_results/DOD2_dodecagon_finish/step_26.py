code = r'''
import mpmath as mp, pickle, time
mp.mp.dps = 620
t0 = time.time()
D1 = '../DOD1_dodecagon/'
r3 = mp.sqrt(3)
X = {1: 1-2*r3/3, 2: 1-r3, 5: 1-r3/3, 6: 1-r3/2}
IDX = [1,2,5,6]
atoms = pickle.load(open(D1+'atoms620.pkl','rb'))
T3a = atoms['T3a']; SwapA = atoms['Swap']
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
    # = Li_{2,1}(b/a, a) (real-analytic); identity: C = -ln(1-a)[Li2(b)-Li2(b/a)] + swap + Li3(b)
    # i.e. swap = C + ln(1-a)[Li2(b)-Li2(b/a)] - Li3(b)  (DOD1 step_45 convention, L1a=-ln(1-a))
    return Creal(a,b) + mp.log(abs(1-a))*(Li2R(b) - Li2R(b/a)) - Li3R(b)
def Li21ser(x,y):
    tot=mp.mpf(0); ps=mp.mpf(0); xp=mp.mpf(1); yp=mp.mpf(1)
    for m in range(1,200000):
        xp*=x
        if m>=2: tot += xp*ps/m**2
        yp*=y; ps += yp/m
        if m>50 and abs(xp*ps)<mp.mpf(10)**(-660): break
    return tot

# CHECK 0: SwpR(x_i,x_j) vs atoms620 Swap (all 12 pairs, includes x-ratio>1 cases)
worst0 = mp.mpf(0)
for (i,j) in sorted(SwapA):
    d = abs(SwpR(X[i],X[j]) - SwapA[(i,j)])
    worst0 = max(worst0,d)
print("CHECK0 SwpR(x_i,x_j) vs atoms620 Swap, worst:", mp.nstr(worst0,4), "t=%ds"%(time.time()-t0))

# CHECK 1: 24 gam/del swaps; series check where convergent
SwG = {}; SwD = {}; worst_ser = mp.mpf(0); nser = 0
for (i,j,k) in sorted(T3a.keys()):
    a1,b1 = gam[(i,j)], gam[(i,k)]
    a2,b2 = dl[(i,j)], dl[(i,k)]
    sg = SwpR(a1,b1); sd = SwpR(a2,b2)
    SwG[(i,j,k)] = sg; SwD[(i,j,k)] = sd
    for (a,b,sw) in [(a1,b1,sg),(a2,b2,sd)]:
        if abs(a)<1 and abs(b/a)<1:
            d = abs(sw - Li21ser(b/a,a)); worst_ser = max(worst_ser,d); nser += 1
print("CHECK1 gam/del swaps vs series (%d convergent cases), worst:"%nser, mp.nstr(worst_ser,4), "t=%ds"%(time.time()-t0))
for k in sorted(SwG): print(k, "SwG=%s SwD=%s"%(mp.nstr(SwG[k],10), mp.nstr(SwD[k],10)))
pickle.dump({'gam':gam,'dl':dl,'SwG':SwG,'SwD':SwD}, open('gdswaps620.pkl','wb'))
print("saved, elapsed", time.time()-t0)
'''
open('step_03.py','w').write(code)
exec(code)
