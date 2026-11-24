# s18: are the submitted certificates SENSITIVE?  Patch one coefficient inside tow3lib and re-run the submitted checkers; they must say False.
import sys, io, contextlib
sys.dont_write_bytecode = True
import sympy as sp
import tow3lib as T3
import prove_V56 as PV
import prove_adj56 as PA
import prove_E3b56 as PE
d, m = T3.d, T3.m
def quiet(f, *a):
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            r = f(*a)
        return ('returned', r, buf.getvalue().strip()[-90:])
    except AssertionError:
        return ('AssertionError', None, buf.getvalue().strip()[-90:])
print('control  non-adjacent (5,0):', quiet(PV.identity, 5, 0))
orig51 = T3._EF[(5,1)]; orig61 = T3._EF[(6,1)]; orig62 = T3._EF[(6,2)]
T3._EF[(5,1)] = orig51 + 1/(12*T3.Pi(2))
print('mutant e(5,1)+const: non-adjacent (5,0):', quiet(PV.identity, 5, 0))
print('mutant e(5,1)+const: non-adjacent (6,1):', quiet(PV.identity, 6, 1))
T3._EF[(5,1)] = orig51
T3._EF[(6,2)] = orig62 + 1/(2*T3.Pi(3))
print('mutant e(6,2)+const: non-adjacent (6,2):', quiet(PV.identity, 6, 2))
T3._EF[(6,2)] = orig62
T3._EF[(6,1)] = orig61 + m**2/(96*T3.Pi(3))
print('mutant e(6,1)+m^2 : non-adjacent (6,0):', quiet(PV.identity, 6, 0))
T3._EF[(6,1)] = orig61
# adjacent + E3 sensitivity (cheap ones)
orig_top = T3.e_top
T3.e_top = lambda delta: orig_top(delta)*sp.Rational(3,2) if delta == 5 else (orig_top(delta)*sp.Rational(4,5) if delta == 6 else orig_top(delta))
T3._a = {} if False else T3._a
print('mutant tops (e(5,2)*3/2, e(6,3)*4/5): adjacent (5,2):', quiet(PA.main, 5, 2))
print('mutant tops: adjacent (6,3):', quiet(PA.main, 6, 3))
T3.e_top = orig_top
print('restored: adjacent (6,3):', quiet(PA.main, 6, 3))
