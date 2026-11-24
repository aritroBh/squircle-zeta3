step_45 = r'''
import pickle, sympy as sp, time
t0=time.time()
E3 = pickle.load(open('E3_fixed.pkl','rb'))
Jtab = pickle.load(open('Jtab_t.pkl','rb'))
A0,A1,B,kap,L2,q0,V02,lA0,lA1,lB = sp.symbols('A0 A1 B kap L2 q0 V02 lA0 lA1 lB')
Jx = sp.Symbol('x')
JL1,JL2,JL3,JAA,Jz2 = sp.symbols('JL1 JL2 JL3 JAA Jz2')
ren = {sp.Symbol('L1'):JL1, sp.Symbol('L2'):JL2, sp.Symbol('L3'):JL3, sp.Symbol('AA'):JAA, sp.Symbol('z2'):Jz2, sp.Symbol('x'):Jx}
# atom coefficient extractors
cA = {}; cL = {}
for key, expr in Jtab.items():
    ex = sp.expand(expr.subs(ren, simultaneous=True))
    cA[key] = ex.coeff(JAA)
    cL[key] = ex.coeff(JL3)
symtab = {s.name: s for s in [A0,A1,B,kap,L2,q0,V02,lA0,lA1,lB]}
print("setup %.1fs"%(time.time()-t0))
PhiA = sp.Integer(0); PhiL = sp.Integer(0)
for (n,a,b,c),co in E3['nl']:
    PhiA += co*cA[(n,a,b,c)]
    PhiL += co*cL[(n,a,b,c)]
print("summed %.1fs"%(time.time()-t0))
PhiA = sp.expand(PhiA); PhiL = sp.expand(PhiL)
print("expanded %.1fs"%(time.time()-t0))
print("PhiA terms:", len(sp.Poly(PhiA, Jx).terms()) if PhiA!=0 else 0)
# factor out sqrt(L2)*B maybe
import pickle as pk
pk.dump({'PhiA':PhiA,'PhiL':PhiL}, open('Phi_syms.pkl','wb'))
print("saved Phi_syms.pkl  %.1fs"%(time.time()-t0))
# print raw sizes
print("PhiA:", PhiA)
print()
print("PhiL:", PhiL)
'''
open('step_45.py','w').write(step_45)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python','step_45.py'],
                   capture_output=True, text=True, timeout=1800)
print(r.stdout)
print("STDERR:", r.stderr[-2000:] if r.stderr else "")
