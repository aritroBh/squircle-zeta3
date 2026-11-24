# A3-mut: SYMBOLIC mutation tests of the closing conditions (the closing test must be able to FAIL).
import sys, time, pickle
sys.dont_write_bytecode=True
W='<project>/research/08_theory/TOW4_tower_depth8'
sys.path.insert(0, W)
import sympy as sp
import core, tow2lib as T2
d, m = T2.d, T2.m; I = sp.I
def tosym(g): return sp.Rational(int(g.r.p),int(g.r.q)) + I*sp.Rational(int(g.i.p),int(g.i.q))
TAU1=[tosym(z) for z in core.TAU1]; TAU2=[tosym(z) for z in core.TAU2]
amine0 = {k: sp.sympify(v, locals={'d':d,'m':m,'I':I}) for k,v in pickle.load(open('a_mine.pkl','rb')).items()}
def run(name, tau1, tau2, amine, maxdelta=8):
    F={}
    def Fget(dp,s): return sp.Integer(0) if s>dp//2 else F[(dp,s)]
    def Lam(delta,s):
        tot=sp.Integer(0)
        for dp in range(delta+1):
            fA=Fget(dp,s).subs({d:d+1,m:m+1},simultaneous=True); fB=Fget(dp,s).subs({d:d+1,m:m},simultaneous=True)
            tot += I*(2*d-1)/(2*d-dp-1)*(tau1[delta-dp]*fA - tau2[delta-dp]*fB)
        return tot
    res={}
    for delta in range(maxdelta+1):
        F[(delta,0)]=sp.cancel(amine[delta])
        for s in range(delta//2+1):
            nxt=sp.cancel(sp.together(Lam(delta,s)-Fget(delta,s)))
            if s+1<=delta//2: F[(delta,s+1)]=nxt
            else: res[delta]= (nxt==0)
    first_fail = [dl for dl in sorted(res) if not res[dl]]
    print(f'{name:60s} closing OK at all delta<=8: {all(res.values())};  failing deltas: {first_fail}', flush=True)
t0=time.time()
run('baseline (my a_delta, my tau)', TAU1, TAU2, amine0)
t1=list(TAU1); t1[8]=t1[8]+sp.Rational(1,128)
run('tau1_8 += 1/128', t1, TAU2, amine0)
t1=list(TAU1); t1[7]=t1[7]+I*sp.Rational(1,64)
run('tau1_7 += i/64', t1, TAU2, amine0)
t2=list(TAU2); t2[8]=-t2[8]
run('tau2_8 -> -tau2_8', TAU1, t2, amine0)
a1=dict(amine0); a1[8]=a1[8]+sp.Rational(1,1000)/(2*d-3)
run('a_8 += 1/(1000(2d-3))', TAU1, TAU2, a1)
a2=dict(amine0); a2[7]=a2[7]+I*sp.Rational(1,1000)*m/((2*d-3)*(2*d-5))
run('a_7 += i m/(1000(2d-3)(2d-5))', TAU1, TAU2, a2)
a3=dict(amine0); a3[8]=a3[8].subs({m:d-m},simultaneous=True)
run('a_8(d,m) -> a_8(d,d-m)', TAU1, TAU2, a3)
a4=dict(amine0); a4[8]=a4[8].subs(I,-I)
run('a_8 -> conj(a_8)', TAU1, TAU2, a4)
print('time', time.time()-t0)
