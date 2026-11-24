
import pickle, numpy as np
from fractions import Fraction as Fr
D = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/solve81.pkl','rb'))
li4v, prodv, Rv = D['li4v'], D['prodv'], D['Rv']
p = (1<<31)-1
def rank_modp(vecs):
    A = [np.array([x % p for x in v], dtype=np.int64) for v in vecs]
    basis = []
    for v in A:
        for b, pw in basis:
            if v[pw]: v = (v - int(v[pw])*b) % p
        nz = np.nonzero(v)[0]
        if len(nz):
            pw = int(nz[0]); basis.append(((v*pow(int(v[pw]), p-2, p)) % p, pw))
    return len(basis)
std8 = ['c', '-c', '1 - c', '1/(c + 1)', 'c/2 + 1/2', '(1 - c)/(c + 1)', '(c - 1)/(c + 1)', '2*c/(c + 1)']
extra = ['c**2', '(1 - c)*(c + 1)', 'c/(c + 1)', '-c**2/((1 - c)*(c + 1))', '1/(1 - c)', '2/(c + 1)', '(c + 1)/c', '-c/(1 - c)', '4*c/(c + 1)**2', '-2*c/(1 - c)']
prodl = list(prodv.values())
for S, label in [(std8, 'std8'), (std8+extra[:4], 'std8+4'), (std8+extra, 'std8+10')]:
    rl = rank_modp([li4v[k] for k in S])
    ra = rank_modp([li4v[k] for k in S] + prodl)
    ok = all(rank_modp([li4v[k] for k in S] + prodl + [Rv[E]]) == ra for E in (-1,0,1))
    print(f"{label}: rank Li4={rl}, rank with products={ra}, R_E in span: {ok}")
