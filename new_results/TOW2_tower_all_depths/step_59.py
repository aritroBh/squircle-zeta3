import sys, os
work = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, work); os.chdir(work)
# append simplified tau table to tow2lib
add = '''

# ---------------- simplified tau constants (all in Q(i)), verified against tau1/tau2 defs ----------------
_tau_tab = {
 0: ((-1-sp.I)/2, (-1+sp.I)/2),
 1: (sp.Rational(3,4)*(1-sp.I), -sp.Rational(3,4)*(1+sp.I)),
 2: (1+sp.Rational(3,4)*sp.I, 1-sp.Rational(3,4)*sp.I),
 3: (sp.Rational(-5,8)+sp.Rational(5,4)*sp.I, sp.Rational(5,8)+sp.Rational(5,4)*sp.I),
 4: (sp.Rational(-23,16)-sp.Rational(7,16)*sp.I, sp.Rational(-23,16)+sp.Rational(7,16)*sp.I),
 5: (sp.Rational(7,32)-sp.Rational(49,32)*sp.I, sp.Rational(-7,32)-sp.Rational(49,32)*sp.I),
 6: (sp.Rational(49,32), sp.Rational(49,32)),
 7: (sp.Rational(3,16)+sp.Rational(93,64)*sp.I, sp.Rational(-3,16)+sp.Rational(93,64)*sp.I),
}
def tau1f(n):
    return _tau_tab[n][0]
def tau2f(n):
    return _tau_tab[n][1]
def check_tau(nmax=7):
    ok = True
    for n in range(0, nmax+1):
        ok &= bool(sp.simplify(_tau_tab[n][0] - tau1(n)) == 0)
        ok &= bool(sp.simplify(_tau_tab[n][1] - tau2(n)) == 0)
    return ok
'''
txt = open('tow2lib.py').read()
open('tow2lib.py','w').write(txt + add)
import importlib, tow2lib
importlib.reload(tow2lib)
print('tau table check:', tow2lib.check_tau())
