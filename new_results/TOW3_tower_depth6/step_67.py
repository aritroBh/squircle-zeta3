import sys, os, json, pickle
from fractions import Fraction as Fr
from math import comb
work = '<project>/research/08_theory/TOW3_tower_depth6'
sys.path.insert(0, work); os.chdir(work)
os.nice(10)
H = pickle.load(open('<project>/research/08_theory/TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb'))
def u(i,j,k,delta):
    return tuple(H[(i,j)][delta][k]) if (i,j) in H and delta in H[(i,j)] and k in H[(i,j)][delta] else (Fr(0),)*4
def qsi_add(x,y): return tuple(a+b for a,b in zip(x,y))
def qsi_scale(x,f): return tuple(a*f for a in x)
def gmul(x,y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def gdiv(x,y):
    den = y[0]**2+y[1]**2
    return ((x[0]*y[0]+x[1]*y[1])/den, (x[1]*y[0]-x[0]*y[1])/den)
def qsi_div(x, y):
    A=(x[0],x[2]); B=(x[1],x[3]); C=(y[0],y[2]); D=(y[1],y[3])
    C2m2D2 = (gmul(C,C)[0]-2*gmul(D,D)[0], gmul(C,C)[1]-2*gmul(D,D)[1])
    numA = (gmul(A,C)[0]-2*gmul(B,D)[0], gmul(A,C)[1]-2*gmul(B,D)[1])
    numB = (gmul(B,C)[0]-gmul(A,D)[0], gmul(B,C)[1]-gmul(A,D)[1])
    return (gdiv(numA, C2m2D2), gdiv(numB, C2m2D2))
def e_ratio(key, delta, s):
    i,j = key
    k0 = (2*i+1) if i==j else i+j
    out = (Fr(0),)*4
    for t in range(0, s+1):
        out = qsi_add(out, qsi_scale(u(i,j,k0+t,delta), Fr((-1)**(s-t)*comb(s,t))))
    return qsi_div(out, u(i,j,k0,0))

# Gaussian+sqrt2 arithmetic: element = (gr, gs) with gr,gs Gaussian rationals; value gr + gs*sqrt2
import sympy as sp
aq = 1+sp.sqrt(2); bq = sp.sqrt(2)-1
ii = sp.I
def wm(mm):
    # (-a-i)^{-m} - (b-i)^{-m} as exact (re + i im) with sqrt2: return sympy
    return (-aq-ii)**(-mm) - (bq-ii)**(-mm)

def beta_minus_a(key, d, M):
    """beta_{d,M} at pole -a from the tower map; returns (rat, s2) Gaussian-rational pair."""
    ikey, j = key
    n = d+1
    parent = (ikey+1, j)
    e00 = u(parent[0], parent[1], n, 0)   # QSI tuple
    E = ((e00[0], e00[2]), (e00[1], e00[3]))  # (gr, gs)
    # Y = sum_m (-1/m) e_ratio(2n-2-(m+1), M-1) w_m ; only m with floor((2d-m-1)/2) >= M-1
    # beta = -(i+1) * 2 Re_i[ E * Y ]
    Yr = (Fr(0),Fr(0)); Ys = (Fr(0),Fr(0))
    for mm in range(1, 2*d):
        deltap = 2*n-2-(mm+1)
        if deltap < 0: break
        if M-1 > deltap//2: continue
        er = e_ratio(parent, deltap, M-1)
        assert er[1] == (Fr(0),Fr(0)), ('sqrt2 in e_ratio', key, deltap, M-1)
        erg = er[0]  # Gaussian rational
        w = sp.expand(wm(mm))
        wr = (sp.Rational(sp.re(w).subs(sp.sqrt(2),0)), sp.Rational(sp.im(w).subs(sp.sqrt(2),0)))
        # split w into rat/s2 parts: w = (wr0 + ws0 s2) + i(wr1 + ws1 s2)
        wre = sp.re(w); wim = sp.im(w)
        wr0 = sp.Rational(wre.subs(sp.sqrt(2),0)); ws0 = sp.Rational(sp.simplify((wre - wr0)/sp.sqrt(2)))
        wr1 = sp.Rational(wim.subs(sp.sqrt(2),0)); ws1 = sp.Rational(sp.simplify((wim - wr1)/sp.sqrt(2)))
        # Y += (-1/m) * erg * w
        f = Fr(-1, mm)
        Yr = (Yr[0] + f*(erg[0]*wr0 - erg[1]*wr1), Yr[1] + f*(erg[0]*wr1 + erg[1]*wr0))
        Ys = (Ys[0] + f*(erg[0]*ws0 - erg[1]*ws1), Ys[1] + f*(erg[0]*ws1 + erg[1]*ws0))
    # E*Y: (Er + Es s2)(Yr + Ys s2) = (ErYr + 2 EsYs) + (ErYs + EsYr) s2
    def ggmul(x,y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
    Zr = (ggmul(E[0],Yr)[0] + 2*ggmul(E[1],Ys)[0], ggmul(E[0],Yr)[1] + 2*ggmul(E[1],Ys)[1])
    Zs = (ggmul(E[0],Ys)[0] + ggmul(E[1],Yr)[0], ggmul(E[0],Ys)[1] + ggmul(E[1],Yr)[1])
    # conj_i(Z) flips im parts
    cZr = (Zr[0], -Zr[1]); cZs = (Zs[0], -Zs[1])
    # 2 Re_i[Z] = Z + conj(Z)
    Rr = (Zr[0]+cZr[0], Zr[1]+cZr[1]); Rs = (Zs[0]+cZs[0], Zs[1]+cZs[1])
    # beta = -(i+1) * (Rr + Rs s2)
    def mul_mi1(x):  # -(i+1)*(x0 + x1 i) = -(x0 - x1) - (x0 + x1) i
        return (-(x[0]-x[1]), -(x[0]+x[1]))
    br = mul_mi1(Rr); bs = mul_mi1(Rs)
    return (br, bs)

tab = json.load(open('<project>/research/08_theory/SUB1_subleading_table/subtop_table.json'))
for d in [2,3,4,5]:
    for key in [(d-j, j) for j in range(0, d//2+1)]:
        top = beta_minus_a(key, d, d)
        sub = beta_minus_a(key, d, d-1)
        exp_top = Fr(tab[str(d)][f'{key}|-a|0']['betas'][-1])
        exp_sub = Fr(tab[str(d)][f'{key}|-a|0']['betas'][-2])
        ok_t = top[0]==(exp_top,Fr(0)) and top[1]==(Fr(0),Fr(0))
        ok_s = sub[0]==(exp_sub,Fr(0)) and sub[1]==(Fr(0),Fr(0))
        print(f'd={d} key={key}: top map={top[0][0]} exp={exp_top} {"OK" if ok_t else "FAIL"}; sub map={sub[0][0]} exp={exp_sub} {"OK" if ok_s else "FAIL"}; s2 parts {top[1]},{sub[1]}')
