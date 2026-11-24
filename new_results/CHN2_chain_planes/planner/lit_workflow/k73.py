import sys, pickle, time; sys.dont_write_bytecode=True
sys.path.insert(0,'.')
import sympy as sp
a1,a2,a3,a4=sp.symbols('a1 a2 a3 a4')
Ms=pickle.load(open('M357_a.pkl','rb')); M3,M5,M7=Ms[3],Ms[5],Ms[7]
pref=a1*a2*a3*a4*(a1+a2+a3+a4)
def K(Mj,Mk): return sp.cancel(sp.expand(Mj[0,1]*Mk[1,0]-Mk[0,1]*Mj[1,0])/pref)
K53=sp.factor(K(M5,M3)); K73=sp.factor(K(M7,M3)); K75=sp.factor(K(M7,M5))
print('K53 deg',sp.Poly(K53,a1,a2,a3,a4).total_degree(),'=',K53)
print('K73 factors:'); 
for f,m in sp.factor_list(K73)[1]: print('   deg',sp.Poly(f,a1,a2,a3,a4).total_degree(),'mult',m, (f if sp.Poly(f,a1,a2,a3,a4).total_degree()<=3 else ''))
print('K75 factors:');
for f,m in sp.factor_list(K75)[1]: print('   deg',sp.Poly(f,a1,a2,a3,a4).total_degree(),'mult',m, (f if sp.Poly(f,a1,a2,a3,a4).total_degree()<=3 else ''))
# does C3 divide K73? K73 mod C3 via polynomial remainder in a1
C3=sp.Poly(K53*(-1) if True else K53,a1)
print('K73 vanishes on Pi1?',sp.simplify(K73.subs({a2:a1,a4:a3},simultaneous=True))==0,' Pi2?',sp.simplify(K73.subs({a3:a1,a4:a2},simultaneous=True))==0)
print('K75 vanishes on Pi1?',sp.simplify(K75.subs({a2:a1,a4:a3},simultaneous=True))==0,' Pi2?',sp.simplify(K75.subs({a3:a1,a4:a2},simultaneous=True))==0)
pickle.dump((K53,K73,K75),open('K.pkl','wb'))
