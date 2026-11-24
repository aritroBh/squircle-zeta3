src = open('fast_cmzv.py').read()
old = """def _aform(j):
    a = -ZET ** j / (1 - ZET ** j)
    return acb_poly([a ** (m + 1) for m in range(N_ + 1)])"""
new = """def _aform(j):
    # omega_j(1-s) = -a_j ds/(1 - a_j s), a_j = -zeta^j/(1-zeta^j): series -a_j^{m+1} s^m
    a = -ZET ** j / (1 - ZET ** j)
    return acb_poly([-a ** (m + 1) for m in range(N_ + 1)])"""
assert old in src
open('fast_cmzv.py','w').write(src.replace(old, new))
import sys; sys.dont_write_bytecode=True; sys.path.insert(0,'.')
import fast_cmzv as FC
FC.setup(62); FC.finish_setup()
import time, mpmath as mp
mp.mp.dps = 65
tests = [((3,), -mp.log(2)), ((1,), mp.pi*mp.j/3), ((-1,0), mp.zeta(2)), ((-1,-1,0), mp.zeta(3)),
         ((-1,-1,-1,0), mp.zeta(4)), ((2,), -mp.log(3)/2 + mp.pi*mp.j/6), ((4,), -mp.log(3)/2 - mp.pi*mp.j/6),
         ((5,), -mp.pi*mp.j/3)]
for w, expect in tests:
    got = FC.cmzv_value(w)
    err = abs(mp.mpc(got.real.mid(), got.imag.mid())-expect)
    print(w, 'err=%.2e'%float(err), 'OK' if err < mp.mpf(10)**(-58) else 'FAIL')
# cross-validate vs mpmath evaluator on random words
import numeric_cmzv as NC, random
NC.set_dps(62)
random.seed(1)
worst=0
for _ in range(12):
    w = tuple(random.choice([-1,0,1,2,3,4,5]) for _ in range(4))
    if w[0]==0 or w[-1]==-1: continue
    a = FC.cmzv_value(w); b = NC.cmzv_value(w)
    err = abs(mp.mpc(a.real.mid(), a.imag.mid())-b)
    worst = max(worst, float(err))
print('cross-check vs mpmath evaluator, worst err over 12 random words: %.2e'%worst)
