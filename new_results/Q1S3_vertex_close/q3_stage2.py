"""Stage 2: build exact Acc expressions for V1 and V3' (y=1/X form).
V1 = -G1 - b*G2 + a*G3
  G1 = int_0^1 [w0 + (r2/2) wa + (-r2/2) w_{-b}] P12        (P12 = g1*g2)
  G2 = a^2 int P12' wa + (1/b) [P12]_0^1                    (a^2 = -(1/b)(-a))
  G3 = b^2 int P12' w_{-b} - (1/a) [P12]_0^1                (b^2 = (1/a)*b)
V3' = -H1 + (b r2/2) H2 + (a r2/2) H3 + r2 b L H4 + r2 a L H5
  H1 = int [w0 - wa - w_{-b}] P13                           (P13 = g1*g3)
  H2 = a^2 int P13' wa + (1/b) [P13]_0^1
  H3 = b^2 int P13' w_{-b} - (1/a) [P13]_0^1
  H4 = g1(1)/(1-b) + sum_{(c,sF) in g1'} a*c int y/((y-sF)(y-a))
  H5 = g1(1)/(1+a) + sum_{(c,sF) in g1'} (-b*c) int y/((y-sF)(y+b))
kernel identities: 1/(1-by) = -a/(y-a); 1/(1+ay) = b/(y+b); y/(1-by) = -a*y/(y-a);
y/(1+ay) = b*y/(y+b).
"""
from q3_engine import *
from q1s_exact import *

r2_ = r2; a_2 = fmul(a_,a_); b_2 = fmul(b_,b_)
hr2 = fscale(r2_, Fr(1,2))          # sqrt2/2

def build_V1():
    acc=Acc()
    P12 = le_prod(g1,g2); dP12 = prod_deriv(P12)
    # G1
    G1=Acc()
    int_simple_kernel(G1, ONE, s0, P12)
    int_simple_kernel(G1, hr2, sa, P12)
    int_simple_kernel(G1, fneg(hr2), smb, P12)
    acc.add_into(G1.scale(fneg(ONE)))
    # G2
    G2=Acc()
    int_formP(G2, a_2, sa, dP12)
    boundary_P_at1(G2, fdiv(ONE,b_), P12)
    acc.add_into(G2.scale(fneg(b_)))
    # G3
    G3=Acc()
    int_formP(G3, b_2, smb, dP12)
    boundary_P_at1(G3, fneg(fdiv(ONE,a_)), P12)
    acc.add_into(G3.scale(a_))
    return acc

def _H45(kernel):
    """kernel 'b': H4; kernel 'a': H5."""
    acc=Acc()
    if kernel=='b':
        den=fsub(ONE,b_); cc0=sa; sk=sa
    else:
        den=fadd(ONE,a_); cc0=fneg(b_); sk=smb
    # boundary g1(1)/den  (note: ('G',0) factor has value 0 at 1)
    for f,c in g1.items():
        if f[0]=='G' and feq(f[1],ZERO): continue
        for lv1 in atom_at1(f):
            acc.add_prod((lv1,), fdiv(c,den))
    for c,sF in le_deriv_terms(g1):
        int_y_rr(acc, fmul(cc0,c), sF, sk)
    return acc

def _append_lv(acc, lv, cF):
    """multiply Acc by F-constant cF and append lv to every product monomial."""
    A=Acc()
    for k,v in acc.P.items():
        A.P[tuple(sorted(k+(lv,),key=repr))]=fmul(v,cF)
    for k,v in acc.W.items():
        w,lvm=k
        A.W[(w,tuple(sorted(lvm+(lv,),key=repr)))]=fmul(v,cF)
    if not feq(acc.R,ZERO):
        A.P[(lv,)]=fadd(A.P.get((lv,),ZERO), fmul(acc.R,cF))
    return A

def build_V3p():
    acc=Acc()
    P13 = le_prod(g1,g3); dP13 = prod_deriv(P13)
    H1=Acc()
    int_simple_kernel(H1, ONE, s0, P13)
    int_simple_kernel(H1, fneg(ONE), sa, P13)
    int_simple_kernel(H1, fneg(ONE), smb, P13)
    acc.add_into(H1.scale(fneg(ONE)))
    H2=Acc()
    int_formP(H2, a_2, sa, dP13)
    boundary_P_at1(H2, fdiv(ONE,b_), P13)
    acc.add_into(H2.scale(fmul(fdiv(b_,TWO),r2_)))   # b*r2/2
    H3=Acc()
    int_formP(H3, b_2, smb, dP13)
    boundary_P_at1(H3, fneg(fdiv(ONE,a_)), P13)
    acc.add_into(H3.scale(fmul(fdiv(a_,TWO),r2_)))   # a*r2/2
    H4=_H45('b')
    acc.add_into(_append_lv(H4, ('lv',a_), fmul(r2_,b_)))
    H5=_H45('a')
    acc.add_into(_append_lv(H5, ('lv',a_), fmul(r2_,a_)))
    return acc

if __name__=="__main__":
    import pickle
    V1=build_V1(); V3=build_V3p()
    print("V1: words %d, prods %d, R %s"%(len(V1.W),len(V1.P),fshow(V1.R)))
    print("V3': words %d, prods %d, R %s"%(len(V3.W),len(V3.P),fshow(V3.R)))
    pickle.dump({'V1':V1,'V3p':V3}, open('q3_accs.pkl','wb'))
