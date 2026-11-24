# s10: which perturbations of the INITIAL data e(delta,0)=a_delta are invisible to the whole T1 / Newton system?
import sys, time
sys.dont_write_bytecode = True
from candlib import *
d, m = T3.d, T3.m; I = sp.I; Pi = T3.Pi
base = base_exprs()
def run(name, exprs):
    ex = dict(base); ex.update(exprs)
    c = Cand(ex)
    tested, bad, per = t1_test(c, [5, 6, 7, 8, 9, 10, 11])
    print(f'{name}: cells {tested}, bad {len(bad)}', flush=True)
run('a_6 + 3/Pi_3(d)       (even homogeneous direction c/Pi_s)', {(6,0): base[(6,0)] + 3/Pi(3)})
run('a_6 + (d-5)/Pi_3(d)   (NOT homogeneous: wrong shape)', {(6,0): base[(6,0)] + (d-5)/Pi(3)})
run('a_6 + m/Pi_3(d)       (m-dependent)', {(6,0): base[(6,0)] + m/Pi(3)})
run('a_5 + 3i(d-4)/Pi_2(d) (odd homogeneous direction c(d-s-2)/Pi_s)', {(5,0): base[(5,0)] + 3*I*(d-4)/Pi(2)})
run('a_5 + 3/Pi_2(d)       (wrong shape)', {(5,0): base[(5,0)] + 3/Pi(2)})
run('a_4 + 3/Pi_2(d)       (even homogeneous, depth 4)', {(4,0): base[(4,0)] + 3/Pi(2)})
run('a_3 + 3i(d-3)/Pi_1(d) (odd homogeneous, depth 3)', {(3,0): base[(3,0)] + 3*I*(d-3)/Pi(1)})
run('a_2 + 3/Pi_1(d)       (even homogeneous, depth 2)', {(2,0): base[(2,0)] + 3/Pi(1)})
