import sys, time, pickle
import mpmath as mp
base = "<project>/research/08_theory/"
sys.path.insert(0, base+"OCT9b_octagon_debug")
from e3build import E3_num_factory
from sympy import Rational as R
t0=time.time()
EDGES = [
 (R(1,4),R(1,5),R(1,5),R(3,2),R(3,5),R(-7,12),[0,R(1,5)],[R(1,4),0]),
 (R(1,9),R(1,4),R(1,4),R(1,3),R(1,2),R(1,9),[R(1,4),0],[R(1,3),R(-1,9)]),
 (R(5,9),R(5,8),R(1,8),R(5,3),R(1,4),R(-5,9),[R(1,8),R(-1,4)],[R(1,3),R(-1,9)]),
 (R(1,8),R(1,5),R(1,5),R(1,4),R(3,5),R(1,24),[0,R(-1,5)],[R(1,8),R(-1,4)])]
mp.mp.dps = 340
Q1 = []
for e, ed in enumerate(EDGES):
    f = E3_num_factory(*ed)
    q = mp.quad(f, [0,1])
    Q1.append(q)
    print(f"edge {e}: Earc [0,1] = {mp.nstr(q, 30)}  ({time.time()-t0:.1f}s)")
pickle.dump([mp.nstr(q, 340) for q in Q1], open("quad330_Q1.pkl","wb"))
print("saved quad330_Q1.pkl")
