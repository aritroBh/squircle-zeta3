
import sys, pickle
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import rul1, itertools

def zbar(S):
    w = rul1.stuple_to_word(S)
    return {w: Fr((-1)**len(S))}
def zeta_word(n): return {(0,)*(n-1)+(1,): Fr(-1)}
def li_half(n):
    out = {}
    for c in itertools.product((0, -1), repeat=n-1):
        w = c + (-1,)
        sgn = Fr((-1)**sum(1 for a in c if a == -1))
        out[w] = out.get(w, Fr(0)) + sgn
    return {k: v for k, v in out.items() if v}
cols9 = pickle.load(open('cols9.pkl','rb'))
colpos = {c:i for i,c in enumerate(cols9)}
LI9 = li_half(9)
bad = [w for w in LI9 if w not in colpos]
print("LI9 words not in cols9:", len(bad))
print(bad[:5])
# check one: (-1,0,0,0,0,-1,-1,-1,-1,-1) length?
w = (-1,0,0,0,0,-1,-1,-1,-1,-1)
print("length of the KeyError word:", len(w))
# is li_half(9) generating some length-10 words?
lens = {}
for w2 in LI9:
    lens[len(w2)] = lens.get(len(w2),0)+1
print(lens)
