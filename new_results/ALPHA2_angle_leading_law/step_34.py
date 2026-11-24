import pickle
from fractions import Fraction as Fr
from math import factorial
def qpart(t): return Fr(*t[0])
def fdiff(seq, n):
    for _ in range(n): seq = [seq[i+1]-seq[i] for i in range(len(seq)-1)]
    return seq
def extract_beta(name, kmax):
    data = pickle.load(open(f'<project>/research/08_theory/ALPHA2_angle_leading_law/kern_{name}.pkl','rb'))
    K = data['kernels']; r = Fr(*data['r']); c = Fr(*data['c'])
    tab = {}
    for d in range(1, 8):
        keys = [(i, j) for i in range(d, (d-1)//2, -1) for j in [d-i]]
        for key in keys:
            seqs = {}
            for k in range(d, kmax+1):
                P = K[k][key]
                for rho, cf in P['real'].items(): seqs.setdefault(rho, []).append(qpart(cf))
            for rho, s in seqs.items():
                if len(s) < d + 3: continue
                top = fdiff(list(s[1:]), d-1)
                assert all(x == top[0] for x in top), (name, d, key, rho)
                beta = top[0]/factorial(d-1)
                if beta != 0: tab[(d, key, rho)] = beta
    return c, r, tab

bs = pickle.load(open('<project>/research/08_theory/ALPHA2_angle_leading_law/beta_squircle.pkl','rb'))
sq = {}
for (d,ij,p),v in bs.items():
    sq[(d,ij,{'0':0,'1':1,'-1':-1,'-r':2}[p])] = v[0]

old = pickle.load(open('<project>/research/08_theory/ALPHA2_angle_leading_law/beta_tables.pkl','rb'))
ALL = {}
for name in ['c12','c35','c13','c45']:
    c = dict(c12=Fr(1,2),c35=Fr(3,5),c13=Fr(1,3),c45=Fr(4,5))[name]
    tab = {(d,ij,Fr(*p)):v for (d,ij,p),v in old[name].items()}
    ALL[name] = (c, tab)
for name in ['c23','c57','c14']:
    c, r, tab = extract_beta(name, 16)
    ALL[name] = (c, tab)

Fd = {}   # Fd[d][c] = value
for name,(c,tab) in ALL.items():
    for d in range(1,8):
        ratios = set()
        for (dd,ij,p),v in tab.items():
            if dd!=d: continue
            marker = 2 if p == -((1-c)/(1+c)) else int(p)
            if (d,ij,marker) in sq and sq[(d,ij,marker)]!=0:
                ratios.add(v/sq[(d,ij,marker)])
        if len(ratios)==1:
            Fd.setdefault(d,{})[c] = ratios.pop()
        elif ratios:
            print("MULTIPLE RATIOS", name, d, ratios)
for d in sorted(Fd):
    print(f"d={d}: ", {str(c): Fd[d][c] for c in sorted(Fd[d])})
pickle.dump(Fd, open('<project>/research/08_theory/ALPHA2_angle_leading_law/Fd_data.pkl','wb'))
