import os, sys, pickle, itertools, random, json
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
import symengine as S
from symengine import Q, ZERO, ONE, P_LET, Q_LET, add, scale, shuffle, obstruction, li_symbol, L2, L3, X, Y
import w4_symbol as W
from fractions import Fraction as F

# --- control A: Li_5 at S-units pairs to 0
R2 = Q(0,1)
sunits = {'b=sqrt2-1': Q(-1,1), '-b': Q(1,-1), '1/2': Q(F(1,2)), '1/sqrt2': R2/2,
          '2-sqrt2': Q(2,-1), 'b^2=3-2sqrt2': Q(3,-2), '1/p': ONE/P_LET, 'kappa1=b': Q(-1,1),
          'sqrt2-1 over 2': (Q(-1,1))/2, '(2-sqrt2)/2=-k0': (Q(2,-1))/2}
okA = True
for name, z in sunits.items():
    obs = obstruction(li_symbol(5, z))
    if obs != [ZERO, ZERO]:
        okA = False; print("Li5 FAIL", name, obs)
print("control A: Li_5 at %d S-units pairs to 0: %s" % (len(sunits), okA))

# --- control B: shuffle products pair to 0 (random words over {0,1} alphabet, weights (1,4),(2,3))
random.seed(11)
okB = True
for _ in range(300):
    k = random.choice([1,2])
    w1 = tuple(random.choice([0,1]) for _ in range(k))
    w2 = tuple(random.choice([0,1]) for _ in range(5-k))
    s1 = {w1: Q(random.randint(-4,4), random.randint(-3,3))}
    s2 = {w2: Q(random.randint(-4,4), random.randint(-3,3))}
    if obstruction(shuffle(s1,s2)) != [ZERO, ZERO]:
        okB = False; break
print("control B: 300 random shuffle products pair to 0:", okB)

# --- control C: M4c W5 pairs to (84, -38 sqrt2): re-pair from symbol.json (planner-style)
MC = "<project>/research/08_theory/M4c_w5_closed_form/"
symj = json.loads(open(MC+"symbol.json").read())['symbol']
Sm4 = {tuple(map(int,w)): Q(F(p), F(q)) for w,(p,q) in symj.items()}
obs_m4 = obstruction(Sm4)
print("control C: M4c W5 obstruction from symbol.json:", obs_m4, " (expect 84, -38*sqrt2)")

# --- independent planner-style re-pairing of my block symbols
raw = pickle.load(open("block_symbols.pkl","rb"))
def parse(s): return {tuple(map(int,w)): Q(F(p), F(q)) for w,(p,q) in s.items()}
Pp = {(0,1): F(1), (1,0): F(-1)}
Qq = {(0,0,1): F(1), (0,1,0): F(-2), (1,0,0): F(1)}
Tt = {(1,0,1): F(2), (1,1,0): F(-1), (0,1,1): F(-1)}
def pair(Sd, A, B):
    tot = Q(0)
    for w, c in Sd.items():
        v = A.get(w[:3], F(0))*B.get(w[3:], F(0)) - A.get(w[2:], F(0))*B.get(w[:2], F(0))
        if v: tot += c*v
    return tot
for k in ['local','single','mixed']:
    Sd = parse(raw['blocks'][k])
    print(f"block {k:7s}: pair QP = {pair(Sd,Qq,Pp)},   pair TP = {pair(Sd,Tt,Pp)}")
Sd = parse(raw['total'])
print("block total  : pair QP =", pair(Sd,Qq,Pp), "   pair TP =", pair(Sd,Tt,Pp))
