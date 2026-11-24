import sys, time, pickle
import mpmath as mp
sys.path.insert(0, ".")
sys.path.insert(0, base+"OCT9b_octagon_debug")
from e3build import quad_Earc, E3_num_factory
from sympy import Rational as R
t0=time.time()
EDGES = [
 (R(1,4),R(1,5),R(1,5),R(3,2),R(3,5),R(-7,12),[0,R(1,5)],[R(1,4),0]),
 (R(1,9),R(1,4),R(1,4),R(1,3),R(1,2),R(1,9),[R(1,4),0],[R(1,3),R(-1,9)]),
 (R(5,9),R(5,8),R(1,8),R(5,3),R(1,4),R(-5,9),[R(1,8),R(-1,4)],[R(1,3),R(-1,9)]),
 (R(1,8),R(1,5),R(1,5),R(1,4),R(3,5),R(1,24),[0,R(-1,5)],[R(1,8),R(-1,4)])]
# first: moderate precision, compare subdivisions & methods
mp.mp.dps = 60
for e, ed in enumerate(EDGES):
    f = E3_num_factory(*ed)
    q1 = mp.quad(f, [0,1])
    q2 = mp.quad(f, [0, R(1,3), R(2,3), 1])
    q3 = mp.quad(f, [0,1], method='gauss-legendre')
    print(f"edge {e}: |[0,1]-[0,1/3,2/3,1]| = {mp.nstr(abs(q1-q2),2)}, |[0,1]-GL| = {mp.nstr(abs(q1-q3),2)}  ({time.time()-t0:.1f}s)")
