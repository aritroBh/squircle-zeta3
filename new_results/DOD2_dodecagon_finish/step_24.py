import mpmath as mp, pickle
mp.mp.dps = 620
D1 = '../DOD1_dodecagon/'
r3 = mp.sqrt(3)
X = {1: 1-2*r3/3, 2: 1-r3, 5: 1-r3/3, 6: 1-r3/2}
atoms = pickle.load(open(D1+'atoms620.pkl','rb'))
Li2R = lambda z: mp.re(mp.polylog(2, mp.mpc(z)))
Li3R = lambda z: mp.re(mp.polylog(3, mp.mpc(z)))
def Creal(a,b):
    return mp.quad(lambda t: mp.log(abs(1-a*t))*mp.log(abs(1-b*t))/t, [0,1])
def SwpR(a,b):
    return Creal(a,b) - mp.log(abs(1-a))*(Li2R(b) - Li2R(b/a)) - Li3R(b)
# compare SwpR(x1,x5) with atoms Swap[(1,5)]
sw = SwpR(X[1], X[5])
print("SwpR(x1,x5) =", mp.nstr(sw, 20))
print("Swap[(1,5)] =", mp.nstr(atoms['Swap'][(1,5)], 20))
print("diff:", mp.nstr(abs(sw - atoms['Swap'][(1,5)]),4))
# Now test the series directly
def Li21ser(x,y):
    tot=mp.mpf(0); ps=mp.mpf(0); xp=mp.mpf(1); yp=mp.mpf(1)
    for m in range(1,200000):
        xp*=x
        if m>=2: tot += xp*ps/m**2
        yp*=y; ps += yp/m
        if m>50 and abs(xp*ps)<mp.mpf(10)**(-660): break
    return tot
print("series(x5/x1, x1) =", mp.nstr(Li21ser(X[5]/X[1], X[1]), 20))
# note: Li21ser(x,y) sums x^n y^m / (n^2 m) over n>m>=1.
# my loop: xp=x^m, ps=sum_{k<=m} y^k/k, term = xp*ps/m^2 -> sum_m x^m (sum_{k<m?} ...) -- CHECK indexing
