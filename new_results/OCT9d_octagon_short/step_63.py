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
Q1 = [mp.mpf(s) for s in pickle.load(open("quad330_Q1.pkl","rb"))]
Q2, Q3 = [], []
for e, ed in enumerate(EDGES):
    f = E3_num_factory(*ed)
    q2 = mp.quad(f, [0, R(1,3), R(2,3), 1])
    q3 = mp.quad(f, [0, R(1,2), 1])
    Q2.append(q2); Q3.append(q3)
    print(f"edge {e}: |Q1-Q2| = {mp.nstr(abs(Q1[e]-q2),2)}, |Q1-Q3| = {mp.nstr(abs(Q1[e]-q3),2)}, |Q2-Q3| = {mp.nstr(abs(q2-q3),2)}  ({time.time()-t0:.1f}s)")
pickle.dump([mp.nstr(q,340) for q in Q2], open("quad330_Q2.pkl","wb"))
pickle.dump([mp.nstr(q,340) for q in Q3], open("quad330_Q3.pkl","wb"))
# combined P'''(1)
LEN = {0: mp.sqrt(41)/20, 1: mp.mpf(5)/36, 2: 5*mp.sqrt(13)/72, 3: mp.sqrt(29)/40}
P3q1 = 12*sum(LEN[e]*Q1[e] for e in range(4))
P3q2 = 12*sum(LEN[e]*Q2[e] for e in range(4))
P3q3 = 12*sum(LEN[e]*Q3[e] for e in range(4))
print("P3 quad [0,1]      =", mp.nstr(P3q1, 70))
print("P3 quad [0,1/3..1] =", mp.nstr(P3q2, 70))
print("P3 quad [0,1/2,1]  =", mp.nstr(P3q3, 70))
print("|P3q1-P3q2| =", mp.nstr(abs(P3q1-P3q2),2), " |P3q1-P3q3| =", mp.nstr(abs(P3q1-P3q3),2))
pickle.dump({'P3q330': mp.nstr(P3q1, 340), 'P3q330_b': mp.nstr(P3q2, 340), 'P3q330_c': mp.nstr(P3q3, 340)}, open("P3_quad330.pkl","wb"))
