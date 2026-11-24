# step_03.py -- D12c: u-substitution form, exact split L_q(p)/8 = 2^{-1/p} + eps*H_q(eps),
# and the sup-norm ruler: L_infty(p) = 8*2^{-1/p} exactly.
import mpmath as mp

def Gq_split(eps, q, dps=80):
    """L_q(1/eps)/8 = 2^{-eps} + eps ∫_0^1 u^{eps-1}(1+u)^{-1-eps}[(1+u^{q(1-eps)})^{1/q}-1] du"""
    mp.mp.dps = dps
    eps, q = mp.mpf(eps), mp.mpf(q)
    g = lambda u: u**(eps-1)*(1+u)**(-1-eps)*((1+u**(q*(1-eps)))**(1/q)-1)
    return 2**(-eps) + eps*mp.quad(g, [0, mp.mpf('0.5'), 1], maxdegree=10)

if __name__ == '__main__':
    import json
    mp.mp.dps = 60
    vals = json.load(open('task1_values.json'))
    for key, q in [("3,2",2),("3,3",3),("3,1.5",1.5),("5,4",4)]:
        p = mp.mpf(key.split(',')[0]); known = mp.mpf(vals[key]['value'])
        print(key, "|split - known| =", mp.nstr(abs(Gq_split(1/p,q)*8-known),3))
    # sup-norm ruler: L_inf(p) = 8*2^{-1/p}; two derivations agree
    # (i) geometric: quarter arc has max=|dx| on [0,1/2], gives 8*2^{-1/p}
    # (ii) q->inf of the u-form: 8*eps∫_0^1 u^{eps-1}(1+u)^{-1-eps}du = 8*2F1(1+eps,eps;1+eps;-1) = 8*2^{-eps}
    for p in [3,5,10]:
        e = 1/mp.mpf(p)
        val = 8*e*mp.quad(lambda u: u**(e-1)*(1+u)**(-1-e), [0,1], maxdegree=10)
        print(f"L_inf({p}): u-form {mp.nstr(val,30)} vs 8*2^-1/p {mp.nstr(8*2**(-e),30)} diff {mp.nstr(abs(val-8*2**(-e)),3)}")
