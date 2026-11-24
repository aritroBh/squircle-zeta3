import mpmath as mp, pickle, sys
sys.path.insert(0,'.')
import engine14 as E
mp.mp.dps=60
Xv={i: E.fev(E.X[i],60) for i in E.IDX}
# G-word evaluator: G(a1..an; 1) = int_0^1 dt/(t-a1) G(a2..an; t), via recursive quad
def Geval(word, z=1):
    if not word: return mp.mpf(1)
    a=word[0]; rest=word[1:]
    # handle singularity at t=a in [0,z]: use the fact that for our points a is real outside (0,1) or at endpoints
    f=lambda t: Geval(rest, t)/(t-a)
    return mp.quad(f, [0,z])
# numeric values of atoms
w={i:(lambda tt,i=i: mp.log(1-Xv[i]*tt)) for i in E.IDX}
u={i:(lambda tt,i=i: 1-Xv[i]*tt) for i in E.IDX}
# 1/x values
r={i: 1/Xv[i] for i in E.IDX}
print("1/x:", {i:mp.nstr(r[i],8) for i in E.IDX})
# T3a(1,2,5) = int w_2 w_5/u_1 dt  vs  -(1/x1) G(1/x2,1/x5,1/x1;1)
val = mp.quad(lambda tt: w[2](tt)*w[5](tt)/u[1](tt),[0,1])
gw = -r[1]*Geval([r[2],r[5],r[1]])
print("T3a(1,2,5):", mp.nstr(val,15), " G-form:", mp.nstr(gw,15), " diff", mp.nstr(abs(val-gw),3))
# T3p(1,2) = int w_2^2/u_1 vs -(1/x1) G(1/x2,1/x2,1/x1;1)
val2 = mp.quad(lambda tt: w[2](tt)**2/u[1](tt),[0,1])
gw2 = -r[1]*Geval([r[2],r[2],r[1]])
print("T3p(1,2):", mp.nstr(val2,15), " G-form:", mp.nstr(gw2,15), " diff", mp.nstr(abs(val2-gw2),3))
# T2(1,2) = int w_2/u_1 vs -(1/x1) G(1/x2,1/x1;1)
val3 = mp.quad(lambda tt: w[2](tt)/u[1](tt),[0,1])
gw3 = -r[1]*Geval([r[2],r[1]])
print("T2(1,2):", mp.nstr(val3,15), " G-form:", mp.nstr(gw3,15), " diff", mp.nstr(abs(val3-gw3),3))
