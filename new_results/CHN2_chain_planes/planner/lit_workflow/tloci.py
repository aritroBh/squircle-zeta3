"""Does the circulant/commuting structure extend off s23=0?  Test T_(23) (2<->3 symmetric kinematics) and T_(14) in the full 5-dim space."""
import sys, pickle; sys.dont_write_bytecode=True
sys.path.insert(0,'.')
import sympy as sp
s1,s2,s3,s4,s5=sp.symbols('s1 s2 s3 s4 s5')
G=pickle.load(open('M35_sym.pkl','rb')); M3,M5=G[3],G[5]
C=sp.expand(M5*M3-M3*M5)
loci={'T_(23): s12=s13,s24=s34  (s4=2s1+s2, s5=s2+2s3)':{s4:2*s1+s2,s5:s2+2*s3},
      'T_(14): s12=s24,s13=s34  (s4=s5=s1+s2+s3)':{s4:s1+s2+s3,s5:s1+s2+s3}}
for name,sub in loci.items():
    print('==',name)
    for nm,X in (('M3',M3),('M5',M5),('[M5,M3]',C)):
        Y=X.subs(sub,simultaneous=True).applyfunc(sp.expand)
        print('  ',nm,'entries (11,12,21,22) factored:',[sp.factor(Y[i,j]) for i in range(2) for j in range(2)] if nm=='[M5,M3]' else '')
        print('     circulant (M11==M22 and M12==M21)?',sp.expand(Y[0,0]-Y[1,1])==0, sp.expand(Y[0,1]-Y[1,0])==0)
