"""Finite weight-two check needed to deduce Delta_(2,2)=0 by coassociativity.

All right weight-two cut factors have endpoints/letters in POINTS. Their
single-log primitive contains only the Li2 arguments enumerated below.
Reduce every one, and prove that the resulting four-element right basis
has an injective full symbol, with pi^2 removed in the right factor.
The final weight-four primitive argument is mathematical (Brown, Prop. 5.2
and Borel's rank formula), stated separately in result.md.
"""
from w4_coproduct import *
from itertools import permutations


def main():
    ban_float()
    points = [ZERO,A,-A,B,-B,ONE]
    args = set()
    for x,y,e in permutations(points,3):
        s = (e-y)/(x-y)
        z = ONE-s if s > ZERO else s
        assert z < ONE
        args.add(z)
    print('Distinct right-factor dilogarithm arguments:',len(args),flush=True)
    for z in sorted(args,key=lambda x:(x.p,x.q)):
        red = DE.reduce_li2(z)
        assert all(c.q == 0 for c in red.d.values())
        assert red.li2_keys() in ([],[('Li2',B)]), (z,red)
        print('Li2',z,'=',red)
    # Right factors: pi^2 is zero; the remaining basis has 4 independent symbols.
    basis = [shuffle(log_symbol(Q(2)),log_symbol(Q(2))),
             shuffle(log_symbol(Q(2)),log_symbol(A)),
             shuffle(log_symbol(A),log_symbol(A)),li_symbol(2,B)]
    words = list(product(range(2),repeat=2))
    mat = flint.fmpq_mat([[flint.fmpq(str(s.get(w,ZERO).p)) for s in basis] for w in words])
    assert all(c.q == 0 for s in basis for c in s.values())
    det = mat.det()
    assert det != 0
    print('Symbol matrix (rows 22,2a,a2,aa; columns lam^2,lam*L,L^2,Li2(b)):')
    print(mat)
    print('Exact determinant:',det)
    print('PROVED: right weight-two symbol is injective on every cut factor.')
    print('Together with Delta_(3,1)=0, coassociativity proves Delta_(2,2)=0.')


if __name__ == '__main__':
    main()
