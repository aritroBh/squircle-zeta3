import sympy as sp
s1,s2,s3,s4,s5,al,be=sp.symbols('s1 s2 s3 s4 s5 al be')
s13=s4-s1-s2; s24=s5-s2-s3
P2=sp.Matrix([[-s3*s4+s1*(s3-s5), s13*s24],[s1*s3,(s1+s2)*(s2+s3)-s4*s5]])
m11=s3*(-s1*(s1+2*s2+s3)+s3*s4+s4**2)+s1*s5*(s1+s5)
m12=-s13*s24*(s1+s2+s3+s4+s5)
m21=s1*s3*(s1+s2+s3-2*(s4+s5))
m22=(s2+s3)*((s1+s2)*(s1+s3)-2*s1*s4)-(2*s1*s3-s4**2+2*s2*(s3+s4))*s5+s4*s5**2
M3=sp.Matrix([[m11,m12],[m21,m22]])
# slice s23=0, a=(a1..a4)=(s24,s12,s34,s13): s=(a2,0,a3,a2+a4,a1+a3)
def at(a1,a2,a3,a4): return {s1:a2,s2:0,s3:a3,s4:a2+a4,s5:a1+a3}
planes={'Pi1 (a1=a2,a3=a4)':at(al,al,be,be),'Pi2 (a1=a3,a2=a4)':at(al,be,al,be),'palindromic (x,y,y,x)':at(al,be,be,al),'generic':at(1,2,3,5)}
for n,sub in planes.items():
    for nm,M in (('P2',P2),('M3',M3)):
        Ms=M.subs(sub).applyfunc(sp.expand)
        print(n,nm,'| offdiag equal:',sp.expand(Ms[0,1]-Ms[1,0])==0,'| diag equal:',sp.expand(Ms[0,0]-Ms[1,1])==0,'| [P2,M3] commute (if M3):', end=' ')
        if nm=='M3':
            P=P2.subs(sub); C=(P*Ms-Ms*P).applyfunc(sp.expand); print(C==sp.zeros(2,2))
        else: print()
