"""Check the (3,1) coaction of W4 against the conjecture, retaining pi^2.

The only weight-three integrals needed have two logarithms. Reuse the
paper's exact antiderivatives and classical Li2/Li3 functional equations.
Numerical identity checks in the imported modules are disabled: this script
proves its equalities by exact algebra, not by a tolerance.
"""
from w4_symbol import *
import polylog_exact as PE
import dilog_exact as DE
import kummer
import flint


def exact_sunit(z):
    assert z != ZERO
    j2, m = logvec(z)
    return int(2*j2), int(m)


# Keep every algebraic operation exact, including selection/sorting of atoms.
PE.sunit_log = DE.sunit_log = kummer.sunit_log = exact_sunit
def akey(a):
    if isinstance(a, str):
        return (0, a, F(0), F(0))
    return (1, a[0], a[1].p, a[1].q)
PE._akey = akey
DE._record = lambda name, z, rhs: DE.CERT.append(f'{name}: Li2({z}) = {rhs}')
P, H, L = PE.P, PE.LAM, PE.LL


def logp(z):
    x, y = logvec(z)
    return H*x + L*y


LT_LOG = ({B: ONE, -A: ONE, ZERO: -ONE}, H*(-2))
MT_LOG = ({A: ONE, -B: ONE, ZERO: -ONE}, H*(-2))


def int_pair(l1, l2, z):
    """Exact int_b^a l1*l2*dv/(v-z); endpoint z=a only with l1=l2=lt."""
    (d1, c1), (d2, c2) = l1, l2
    ans = P()
    for lo, hi in [(B, ONE), (ONE, A)]:
        if z in (lo, hi):
            assert z == A and l1 == l2 == LT_LOG
            # lt(a)=0: normalize each logarithm at a before polarizing.
            assert not (c1 + sum((logp(A-x)*c for x,c in d1.items()), P())).d
            for x, cx in d1.items():
                for y, cy in d2.items():
                    ans += PE.int_reg2(x, y, A, lo, hi)*cx*cy
            continue
        for x, cx in d1.items():
            for y, cy in d2.items():
                ans += PE.int_loglog(x, y, z, lo, hi)*cx*cy
        for x, cx in d1.items():
            ans += PE.int_log(x, z, lo, hi)*cx*c2
        for y, cy in d2.items():
            ans += PE.int_log(y, z, lo, hi)*cy*c1
        ans += PE.int_one(z, lo, hi)*c1*c2
    return ans


@cache
def pair(z, mixed=False):
    if z == B:
        if not mixed:
            # FP lt(b)=K. Integrate d(lt^3)/3 and remove the endpoint pole.
            K = L - H*F(1,2)
            return -(K*K*K)*F(1,3) - pair(-A) + pair(ZERO)
        # Integrate d(lt^2*mt)/2; both endpoint values are zero.
        return -pair(-A, True) + pair(ZERO, True) - (
            pair(A) + pair(-B) - pair(ZERO))*F(1,2)
    return int_pair(LT_LOG, MT_LOG if mixed else LT_LOG, z)


def tensor(p, ell):
    return [p*x for x in logvec(ell)]


def vadd(*vs):
    return [sum((v[i] for v in vs), P()) for i in range(2)]


def vscale(v, c):
    return [x*c for x in v]


def lhs_coproduct():
    # Summed delete-one-letter formula, with the inner shuffle telescoped.
    # For int prod(l_j) omega_z:
    # sum_j I_z(prod_except_j) tensor (k_j+sum_x c_jx log|x-z|)
    # -sum_jx c_jx I_x(prod_except_j) tensor log|x-z|
    # +[prod(l_j(e)) tensor log|e-z|]_b^a.
    def pure(z):
        out = [P(), P()]
        d, c = LT_LOG
        C = c + sum((logp(x-z)*cx for x,cx in d.items()), P())
        cv = [C.d.get(('lam',), ZERO), C.d.get(('L',), ZERO)]
        out = vadd(out, [pair(z)*x*3 for x in cv])
        for x, cx in d.items():
            out = vadd(out, vscale(tensor(pair(x), x-z), -3*cx))
        K = L-H*F(1,2)
        return vadd(out, vscale(tensor(K*K*K, B-z), -1))
    def mixed():
        out = [P(), P()]
        for factors, mixed_flag, multiplicity in [(LT_LOG, True, 2), (MT_LOG, False, 1)]:
            d, c = factors
            C = c + sum((logp(x)*cx for x,cx in d.items()), P())
            cv = [C.d.get(('lam',), ZERO), C.d.get(('L',), ZERO)]
            out = vadd(out, [pair(ZERO, mixed_flag)*x*multiplicity for x in cv])
            for x, cx in d.items():
                out = vadd(out, vscale(tensor(pair(x, mixed_flag), x), -multiplicity*cx))
        return out
    return vadd(vscale(pure(ZERO), (8-7*R)/48),
                vscale(pure(-A), Q(F(-1,3))), vscale(mixed(), R/16))


def atom_p(a):
    if a[0] == 'log':
        return logp(a[1])
    if a == ('pi2',):
        return PE.PI2
    if a == ('z3',):
        return PE.Z3
    assert a[0] == 'li'
    n, z = a[1:]
    assert n in (1,2,3)
    return -logp(ONE-z) if n == 1 else (PE.LI2(z) if n == 2 else PE.LI3(z))


def rhs_coproduct():
    out = [P(), P()]
    for c, atoms in conjecture_terms():
        for i, a in enumerate(atoms):
            if a[0] == 'li':
                n, z = a[1:]
                da = tensor(atom_p(('li', n-1, z)), z)
            elif a[0] == 'log':
                da = [P.const(x) for x in logvec(a[1])]
            else:
                assert a in [('pi2',), ('z3',)]
                continue  # pi^2 and zeta(3) have no (w-1,1) component
            for j, b in enumerate(atoms):
                if j != i:
                    da = [v*atom_p(b) for v in da]
            out = vadd(out, vscale(da, c))
    return out


def li3_reduce(expr):
    """Exact Gaussian certificate using only the paper's four Kummer instances."""
    basis = {B, -B, R/2}
    extra = [(f'Kummer x={x}, y={y}', kummer.kummer_relation(x,y))
             for x,y in [(B,R/2), (2-R,R/2), (R/2,2*B), (-ONE,B*B)]]
    args = set(expr.atoms('Li3')) | basis
    for _, rel in extra:
        args |= rel.atoms('Li3')
    rels = {}
    frontier = args.copy()
    # ponytail: four closure rounds suffice here; enlarge only for a new failed certificate.
    for _ in range(4):
        new = set()
        for z in sorted(frontier, key=lambda x: (x.p,x.q)):
            for name, rel in PE.li3_relations(z):
                try:
                    if any(abs(exact_sunit(a)[0])>12 or abs(exact_sunit(a)[1])>8
                           for a in rel.atoms('Li3')):
                        continue
                except (ValueError, AssertionError):
                    continue
                key = repr(rel)
                rels[key] = (f'{name} [seed={z}]', rel)
                new |= rel.atoms('Li3') - args
        args |= new
        frontier = new
        if not frontier:
            break
    allrels = list(rels.values()) + extra + [
        (f'Li3 special value {z}', PE.LI3(z)-val) for z,val in PE.LI3_VALUES.items()]
    assert all(c.q == 0 for _,rel in allrels for c in rel.d.values())
    atoms = sorted(set.union(expr.atoms('Li3'), *(r.atoms('Li3') for _,r in allrels))-basis,
                   key=lambda x: (x.p,x.q))
    rows = []
    for a in atoms:
        row = []
        for _, rel in allrels:
            c = rel.d.get((('Li3',a),), ZERO)
            assert c.q == 0
            row.append(flint.fmpq(str(c.p)))
        c = expr.d.get((('Li3',a),), ZERO)
        rows.append(row + [flint.fmpq(str(c.p)), flint.fmpq(str(c.q))])
    mat, rank = flint.fmpq_mat(rows).rref()
    mu = [ZERO]*len(allrels)
    for i in range(rank):
        pivot = next(j for j in range(mat.ncols()) if mat[i,j])
        assert pivot < len(allrels), 'Li3 relation span insufficient'
        mu[pivot] = Q(F(str(mat[i,len(allrels)])), F(str(mat[i,len(allrels)+1])))
    red = expr
    for c, (name, rel) in zip(mu, allrels):
        if c != ZERO:
            print('  relation multiplier:', c, name)
            print('    exact zero relation:', rel)
            red -= rel*c
    assert not red.atoms('Li3')-basis
    return red


def coaction_shuffle_check():
    """Compare the telescoped (3,1) formula with every literal GPL cut."""
    for z,factors in [(ZERO,(LT,LT,LT)),(-A,(LT,LT,LT)),(ZERO,(LT,LT,MT))]:
        direct, collapsed = {}, {}
        def put(out,e,w,k,ell,cf):
            for i,x in enumerate(ell):
                key = (e if w else ONE,w,k,i)
                out[key] = out.get(key,ZERO)+cf*x
        for (inner,k),cf in log_product(*factors).items():
            w = (z,)+inner
            for e,sgn in [(A,1),(B,-1)]:
                points = (e,)+w+(ONE,)
                for j,x in enumerate(w,1):
                    ell = tuple(a-b for a,b in zip(logvec(x-points[j-1]),logvec(x-points[j+1])))
                    put(direct,e,w[:j-1]+w[j:],k,ell,cf*sgn)
                if k:
                    put(direct,e,w,k-1,(F(1),F(0)),cf*k*sgn)
        for j,f in enumerate(factors):
            lf = LT_LOG if f == LT else MT_LOG
            d,c = lf
            ellp = c+sum((logp(x-z)*cx for x,cx in d.items()),P())
            ell = [ellp.d.get((name,),ZERO) for name in ('lam','L')]
            rest = factors[:j]+factors[j+1:]
            for (w,k),cf in log_product(*rest).items():
                for e,sgn in [(A,1),(B,-1)]:
                    put(collapsed,e,(z,)+w,k,ell,cf*sgn)
                    for x,cx in d.items():
                        put(collapsed,e,(x,)+w,k,logvec(x-z),-cf*cx*sgn)
        for (w,k),cf in log_product(*factors).items():
            for e,sgn in [(A,1),(B,-1)]:
                put(collapsed,e,w,k,logvec(e-z),cf*sgn)
        assert clean(direct) == clean(collapsed), (z, factors)


def main():
    ban_float()
    checks = PE.verify_antiderivatives()
    assert checks == [0]*5
    print('Five exact antiderivative differentiation checks: PASS', flush=True)
    coaction_shuffle_check()
    print('Telescoped (3,1) formula versus literal GPL cuts: PASS',flush=True)
    lhs, rhs = lhs_coproduct(), rhs_coproduct()
    for i, name in enumerate(['log2','log a']):
        diff = PE.reduce_li2_all(lhs[i]-rhs[i])
        print('Reducing coefficient of', name, flush=True)
        red = li3_reduce(diff)
        print('Delta_(3,1)(W4-RHS), coefficient of', name, '=', red, flush=True)
        assert not red.d
    print('PROVED-MODULO-COACTION-KERNEL: Delta_(3,1)(W4-RHS) = 0')


if __name__ == '__main__':
    main()
