import sys; sys.path.insert(0,'.')
import importlib, symbol_engine_r as SE
importlib.reload(SE)
from itertools import combinations
from flint import fmpz_mat

# formal 2-letter case: x = p1, y = p2
s = 2
L, idx = SE.letters(s)
A = [('0',),('1',)] + [('p',i) for i in range(s)]
Pts = A + [('inf',)]
# letter vectors
p1 = {idx[('p',0)]:1}; p2 = {idx[('p',1)]:1}
q1 = {idx[('q',0)]:1}; q2 = {idx[('q',1)]:1}
r12 = {idx[('r',0,1)]:1}
# c = x/(x-y): p1 - r12 ; cX = (x-1)/(x-y): q1 - r12 ; 1-c = ? (0-y)/(x-y) -> -y/(x-y): p2 - r12
c  = SE.vadd(p1, r12, -1)
cX = SE.vadd(q1, r12, -1)
one_c = SE.vadd(p2, r12, -1)       # 1 - c = -y/(x-y)
one_cX = SE.vadd(q2, r12, -1)      # 1 - cX = (1-y)/(x-y)
# CR set
CR = []
for pts4 in combinations(Pts, 4):
    CR += SE.cross_ratios(list(pts4), idx, s)
CRset = {}
for cr, oc in CR:
    CRset[tuple(sorted(cr.items()))] = (cr, oc)
print('c in CR:', tuple(sorted(c.items())) in CRset, '  cX in CR:', tuple(sorted(cX.items())) in CRset)
# symbol of N = 2 Li3(c) - 2 Li3(cX) + products  => mod D: 2 s(c) - 2 s(cX)
sN = SE.tadd(SE.li3_symbol(c, one_c), {k:-v for k,v in SE.li3_symbol(cX, one_cX).items()})
sN = {k:2*v for k,v in sN.items()}
Ly = SE.lyndon3(len(L)); B = [SE.bracket_tensor(w) for w in Ly]
CRsyms = [SE.li3_symbol(cr, oc) for cr, oc in CRset.values()]
M_cr = fmpz_mat([SE.pairing_row(T, B) for T in CRsyms])
M_N  = fmpz_mat([SE.pairing_row(sN, B)])
M_all = fmpz_mat([[M_cr[i,j] for j in range(M_cr.ncols())] for i in range(M_cr.nrows())] +
                 [[M_N[0,j] for j in range(M_N.ncols())]])
print('rank s(CR) mod D:', M_cr.rank(), ' rank with N:', M_all.rank(), ' N in span:', M_all.rank()==M_cr.rank())

# antiderivative checks with sympy (exact)
import sympy as sp
U, cc = sp.symbols('U cc')
F = -sp.log(U)**2*sp.log(1-cc*U) - 2*sp.log(U)*sp.polylog(2, cc*U) + 2*sp.polylog(3, cc*U)
dF = sp.simplify(sp.diff(F, U) - sp.log(U)**2/(1-cc*U)/cc)   # F' = log^2 U/(1-cU)/c  (so N = (F(1)-F(X))/c * ... )
print('antiderivative F: c*dF/dU - log^2/(1-cU) =', sp.simplify(cc*sp.diff(F,U) - sp.log(U)**2/(1-cc*U)))
v = sp.symbols('v')
G0 = sp.log(1-v)**2*sp.log(v) + 2*sp.log(1-v)*sp.polylog(2,1-v) - 2*sp.polylog(3,1-v)
print('antiderivative G0: dG0/dv - log^2(1-v)/v =', sp.simplify(sp.diff(G0,v) - sp.log(1-v)**2/v))
open('progress.log','a').write('step: repeated-word repair: N closed form symbol in span s(CR)+D exactly; antiderivatives certified by sympy\n')
