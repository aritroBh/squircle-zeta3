"""M3,M5,M7 on s23=0 slice (our kinematics) symbolic in a=(a1..a4); commutators on Pi1, Pi2; swap 2<->3 action on a."""
import sys, pickle, time; sys.dont_write_bytecode=True
sys.path.insert(0,'.')
import sympy as sp
from fractions import Fraction as F
import phi
a1,a2,a3,a4=sp.symbols('a1 a2 a3 a4')
def kin_a(a):
    A1,A2,A3,A4=a
    return dict(s12=A3,s13=A1,s23=sp.Integer(0),s24=A4,s34=A2,s123=A1+A3,s234=A2+A4)
def cols_sym(a,Wmax):
    k=kin_a(a); e0,e1=phi.e_mats(**k)
    e0=[[sp.expand(x) for x in r] for r in e0]; e1=[[sp.expand(x) for x in r] for r in e1]
    e={'0':e0,'1':e1}
    cur={'':[sp.Integer(1)]+[sp.Integer(0)]*5}; res={}
    for W in range(1,Wmax+1):
        nxt={}
        for w,v in cur.items():
            for ch in '01':
                M=e[ch]; nxt[w+ch]=[sp.expand(sum(M[i][j]*v[j] for j in range(6))) for i in range(6)]
        cur=nxt
        if W in (3,5,7):
            tot=[0,0]
            for w,v in cur.items():
                if v[0]==0 and v[1]==0: continue
                zv=phi.zeta_vec(w)[0]          # coefficient of single zeta_W (basis index 0)
                if zv!=0:
                    c=sp.Rational(zv.numerator,zv.denominator)
                    tot[0]+=v[0]*c; tot[1]+=v[1]*c
            res[W]=[sp.expand(tot[0]),sp.expand(tot[1])]
    return res
def Fk(a,k):   # F2 linear term at s23=0: s12^k + 0 - s12^k = 0  -> no contribution
    return 0
def Mmats(Wmax=7):
    t=time.time()
    r=cols_sym((a1,a2,a3,a4),Wmax)
    ap=(a3,a4,a1,a2)             # 2<->3 swap acts as a -> (a3,a4,a1,a2) on our s23=0 slice
    rp=cols_sym(ap,Wmax) if False else None
    out={}
    sw={a1:a3,a2:a4,a3:a1,a4:a2}
    for k in (3,5,7):
        A,B=r[k]
        out[k]=sp.Matrix([[A,B],[B.subs(sw,simultaneous=True),A.subs(sw,simultaneous=True)]])
    print('M3,M5,M7 built in %.1fs'%(time.time()-t))
    return out
if __name__=='__main__':
    Ms=Mmats(7)
    pickle.dump(Ms,open('M357_a.pkl','wb'))
    M3,M5,M7=Ms[3],Ms[5],Ms[7]
    print('M3=',[[sp.factor(x) for x in M3.row(i)] for i in range(2)])
    # check M3,M5 here against symbolic general ones restricted
    G=pickle.load(open('M35_sym.pkl','rb'))
    s1,s2,s3,s4,s5=sp.symbols('s1 s2 s3 s4 s5')
    sub={s1:a3,s2:0,s3:a2,s4:a1+a3,s5:a2+a4}
    print('consistency M3:',sp.expand(G[3].subs(sub,simultaneous=True)-M3)==sp.zeros(2,2),' M5:',sp.expand(G[5].subs(sub,simultaneous=True)-M5)==sp.zeros(2,2))
    C53=sp.expand(M5*M3-M3*M5)
    print('[M5,M3]_11 factor:',sp.factor(C53[0,0]))
    for name,sub1 in (('Pi1 (a1=a2,a3=a4)',{a2:a1,a4:a3}),('Pi2 (a1=a3,a2=a4)',{a3:a1,a4:a2})):
        print('=====',name)
        N3=M3.subs(sub1,simultaneous=True); N5=M5.subs(sub1,simultaneous=True); N7=M7.subs(sub1,simultaneous=True)
        for nm,N in (('M3',N3),('M5',N5),('M7',N7)):
            print(nm,'=',[[sp.factor(x) for x in N.row(i)] for i in range(2)])
        for nm,(X,Y) in (('[M5,M3]',(N5,N3)),('[M7,M3]',(N7,N3)),('[M7,M5]',(N7,N5))):
            C=sp.expand(X*Y-Y*X)
            print(nm,'=',[[sp.factor(x) for x in C.row(i)] for i in range(2)],' det=',sp.factor(C.det()))
