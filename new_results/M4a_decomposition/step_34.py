import sys, json, pickle
sys.dont_write_bytecode = True
sys.path.insert(0,'.')
import mpmath as mp
from fractions import Fraction as Fr

mp.mp.dps = 160
r2 = mp.sqrt(2); aN = 1 + r2; bN = r2 - 1
omega = json.load(open('omega4.json'))
terms = {tuple(map(int, k.split(','))): [(z, mp.mpf(c[0]) + r2*mp.mpf(c[1])) for z, c in v]
         for k, v in omega.items()}
pole = {'0': mp.mpf(0), 'a': aN, '-a': -aN, 'b': bN, '-b': -bN}

def piece(i, j):
    ts = terms[(i, j)]
    def f(v):
        l1 = mp.log((v - bN)*(v + aN)/(4*v))
        l2 = mp.log((aN - v)*(v + bN)/(4*v))
        return l1**i * l2**j * mp.fsum(c/(v - pole[z]) for z, c in ts)
    val, err = mp.quad(f, [bN, mp.mpf(1)/2, 1, 2, aN], error=True, maxdegree=14)
    return val, err

for ij in [(1,1),(2,0)]:
    val, err = piece(*ij)
    print(ij, mp.nstr(val, 80), 'err', mp.nstr(err, 3))
W34 = mp.fsum([piece(1,1)[0], piece(2,0)[0]])
tails = pickle.load(open('m4_tails_exact.pkl','rb'))
W3exact = tails['W3p'].num()
print('W34 quad  =', mp.nstr(W34, 110))
print('W3exact   =', mp.nstr(W3exact, 110))
print('diff =', mp.nstr(W34 - W3exact, 5))
