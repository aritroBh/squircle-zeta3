import sys, os, pickle
sys.dont_write_bytecode=True
sys.path.insert(0, os.getcwd())
import chk_forced as C
import tow4lib as T4, sympy as sp
from chk_forced import G, Fr
bad=0;tot=0
for (dl,s) in sorted(pickle.load(open('forced_F78.pkl','rb')).keys()):
    expr = T4.forced(dl,s)   # cached sympy expression in d,m
    for (dd,mm) in [(6,6),(6,4),(7,3),(9,5),(11,2),(13,13)]:
        if 2*dd-dl-1==0: continue
        x = C.ev(sp.sympify(expr), {'d':G(dd),'m':G(mm)})
        y = C.F(dl,s,dd,mm)
        tot+=1
        if not (x==y): bad+=1; print('cache!=own',dl,s,dd,mm)
print('cached forced_F78.pkl vs my own recursion: entries x points =',tot,'bad',bad)
