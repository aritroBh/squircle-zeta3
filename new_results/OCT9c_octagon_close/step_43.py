import sys, time, pickle
sys.path.insert(0,'.')
from oct9c_common import *
D = load_all()
mp.mp.dps = 370
def Li21val(a, c):
    a, c = mp.mpf(a), mp.mpf(c)
    if abs(a) <= 1 and abs(a*c) < 1:
        N = int((mp.mp.dps+20)/(-mp.log10(max(abs(a), abs(a*c))))) + 50
        return mp.fsum([a**n * c**m / (mp.mpf(n)**2 * m) for n in range(2, N) for m in range(1, n)])
    F = mp.quad(lambda t: mp.log(t)*mp.log(1-a*c*t)/(1-a*t), [0,1])
    return a*F

PTS = {0:[(sp.Rational(-3,2),sp.Rational(-7,18)),(sp.Rational(7,12),sp.Rational(-18,7)),
          (sp.Rational(3,5),sp.Rational(-7,3)),(sp.Rational(-7,5),sp.Rational(-3,7))],
       1:[(sp.Rational(-1,3),sp.Rational(1,3)),(sp.Rational(-1,9),sp.Rational(3)),
          (sp.Rational(1,4),sp.Rational(2,5)),(sp.Rational(1,10),sp.Rational(5,2))],
       2:[(sp.Rational(-5,3),sp.Rational(-1,3)),(sp.Rational(5,9),sp.Rational(-3)),
          (sp.Rational(5,8),sp.Rational(-2)),(sp.Rational(-5,4),sp.Rational(-1,2))],
       3:[(sp.Rational(-1,4),sp.Rational(1,6)),(sp.Rational(-1,24),sp.Rational(6)),
          (sp.Rational(1,5),sp.Rational(1,5)),(sp.Rational(1,25),sp.Rational(5))]}
t0=time.time()
mine = {}
for e in range(4):
    for j,(a,c) in enumerate(PTS[e]):
        v = Li21val(a,c)
        mine[f"e{e}L{j+1}"] = mp.nstr(v, 360)
        print(f"e{e}L{j+1} = Li21({a},{c}) done [{time.time()-t0:.0f}s]", flush=True)
pickle.dump(mine, open('lam_mine.pkl','wb'))
print("saved lam_mine.pkl")
