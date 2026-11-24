lib = open('<project>/research/08_theory/SRF2_cube_residue/srf2_lib.py').read()
addition = '''

def J2_int(eps, h=arb(1)/32, T=4.6):
    """J(eps) with TWO subtractions:
       J = 1/(2e+1) - (2+2e)/(2e+2) + int_0^1 a^{2e}[(1+a)^{-2-2e} S0 - 1 + (2+2e) a] da.
       Integrand O(a^{2e+2}), smooth through eps=-1."""
    nodes = ts_nodes(h, T)
    J2 = arb(0)
    for (a, am, wa, pa) in nodes:
        la = a.log(); l1a = (1+a).log()
        A2 = a*a*(arb(-2*eps)*la).exp()
        S0 = (1+A2).sqrt()
        f = S0*expm1_s(arb(-(2+2*eps))*l1a) + A2/(S0+1) + arb(2+2*eps)*a
        J2 += (arb(2*eps)*la).exp()*f*wa
    return 1/(2*eps+1) - arb(2+2*eps)/(2*eps+2) + J2

def W2_int(eps, h=arb(1)/32, T=4.6):
    """W(eps) with the corner monomial removed:
       W = (q(q+1)/2)/(e+2) * [1/(2e+2) - (2+q)/(2e+3) + K2] + int int w * ht,
       ht = h_W - (q(q+1)/2) a^2 c^2 (1+a)^{-2-q} S0,
       K2 = int_0^1 a^{2e+1}[(1+a)^{-2-q} S0 - 1 + (2+q)a] da."""
    q = 1 + 2*eps
    nodes = ts_nodes(h, T)
    K2 = arb(0)
    for (a, am, wa, pa) in nodes:
        la = a.log(); l1a = (1+a).log()
        A2 = a*a*(arb(-2*eps)*la).exp()
        S0 = (1+A2).sqrt()
        f = S0*expm1_s(arb(-(2+q))*l1a) + A2/(S0+1) + arb(2+q)*a
        K2 += (arb(2*eps+1)*la).exp()*f*wa
    K = 1/(2*eps+2) - arb(2+q)/(2*eps+3) + K2
    cpre = []
    for (c, cm, wc, pc) in nodes:
        lc = c.log()
        C2 = c*c*(arb(-2*eps)*lc).exp()
        prefc = (arb(eps)*lc).exp()/c
        cpre.append((c, C2, prefc, wc))
    tot = arb(0)
    for (a, am, wa, pa) in nodes:
        la = a.log(); l1a = (1+a).log()
        A2 = a*a*(arb(-2*eps)*la).exp()
        S0 = (1+A2).sqrt()
        prefa = (arb(2*eps)*la).exp()/a
        outfact = (arb(-2*eps)*l1a).exp()/(1+a)
        corner_a = (arb(q*(q+1)/2))*a*a*(arb(-(2+q))*l1a).exp()*S0
        tot_row = arb(0)
        for (c, C2, prefc, wc) in cpre:
            u = a*c/(1+a)
            S = (1+A2+A2*C2).sqrt()
            br = Phi(u, q)*S + (1 - q*u)*(A2*C2)/(S+S0)
            ht = outfact*br - corner_a*c*c
            tot_row += prefc*ht*wc
        tot += prefa*tot_row*wa
    return tot + (q*(q+1)/2)/(eps+2)*K

def V_cont2(eps, h=arb(1)/32, T=4.6):
    q = 1 + 2*eps
    return W2_int(eps, h, T) - q/(eps+1)*J2_int(eps, h, T)
'''
lib += addition
open('<project>/research/08_theory/SRF2_cube_residue/srf2_lib.py','w').write(lib)
print("extended")
