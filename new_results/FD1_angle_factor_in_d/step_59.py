import sys, pickle; fd1='<project>/research/08_theory/FD1_angle_factor_in_d'
sys.path.insert(0,fd1)
import importlib, fd1_lib; importlib.reload(fd1_lib)
from fd1_lib import squircle_beta
from fractions import Fraction as Fr
from math import factorial, lcm, gcd
def qpart(t): return Fr(*t[0])
def fdiff(seq, n):
    for _ in range(n): seq = [seq[i+1]-seq[i] for i in range(len(seq)-1)]
    return seq
def extract_beta20(name, dmax=10):
    data = pickle.load(open(f'{fd1}/kern20_{name}.pkl','rb'))
    K = data['kernels']; r = Fr(*data['r']); c = Fr(*data['c'])
    tab = {}
    for d in range(1, dmax+1):
        keys = [(i, j) for i in range(d, (d-1)//2, -1) for j in [d-i]]
        for key in keys:
            seqs = {}
            for k in range(d, 21):
                if k not in K or key not in K[k]: continue
                for rho, cf in K[k][key]['real'].items(): seqs.setdefault(rho, []).append(qpart(cf))
            for rho, s in seqs.items():
                if len(s) < d + 1: continue
                top = fdiff(list(s[1:]), d-1)
                if d <= 9: assert all(x == top[0] for x in top), (name, d, key, rho)
                beta = top[0]/factorial(d-1)
                if beta != 0: tab[(d, key, Fr(*rho))] = beta
    return c, r, tab
polmap={-1:'-a',1:'a',0:'0'}
Fd10={}
for name in ['c12','c35','c13','c45','c23','c57','c14','c67']:
    c, r, tab = extract_beta20(name)
    for d in [10]:
        ratios = {}
        for (dd,ij,p),v in tab.items():
            if dd!=d: continue
            pole = '-b' if p == -r else polmap.get(int(p))
            if pole is None: continue
            sb = squircle_beta(d, ij[0], ij[1], pole)
            if sb != 0: ratios[(ij,pole)] = v/sb
        uniq = set(ratios.values())
        assert len(uniq)==1, (name, uniq)
        Fd10.setdefault(d,{})[c]=uniq.pop()
pickle.dump(Fd10, open(fd1+'/Fd_data_d10.pkl','wb'))
F10=Fd10[10]
pts=[((1-c*c)/(c*c), v) for c,v in F10.items()]
pts.append((Fr(1),Fr(1)))
def solve_poly(pts, deg):
    base=pts[:deg+1]
    M=[[x**k for k in range(deg,-1,-1)]+[y] for x,y in base]
    n=deg+1
    for col in range(n):
        piv=next(r for r in range(col,n) if M[r][col]!=0)
        M[col],M[piv]=M[piv],M[col]
        pv=M[col][col]; M[col]=[v/pv for v in M[col]]
        for r in range(n):
            if r!=col and M[r][col]!=0:
                f=M[r][col]; M[r]=[M[r][i]-f*M[col][i] for i in range(n+1)]
    coeffs=[M[r][n] for r in range(n)]
    for x,y in pts:
        assert sum(coeffs[i]*x**(deg-i) for i in range(deg+1))==y,(x,y)
    return coeffs
c10=solve_poly(pts,5)
L=1
for cc in c10: L=lcm(L,cc.denominator)
ints=[int(cc*L) for cc in c10]
g=0
for v in ints: g=gcd(g,abs(v))
ints=[v//g for v in ints]
print("f_10 integer poly (t^5..t^0):", ints, " f_10(1)=",sum(ints))
print("f_10(1) predicted A_10*2^(s2(10)-1) =", 50521*2)
print("LC predicted 10!/2^(10-s2(10)) =", 3628800//2**8)
open(fd1+'/progress.log','a').write('step_56: F_10 re-extracted with corrected squircle table; exact integer poly fit on 9 pts\n')
