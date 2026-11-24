import sys, time
import sympy as sp
sys.path.insert(0, '<project>/research/08_theory/TOW4_tower_depth8')
import tow4lib as T4, tow2lib as T2
d, m = T4.d, T4.m
i = sp.I
tau1 = T4.tau1f; tau2 = T4.tau2f

# descending family ed: ed(delta,top)=zigzag top; ed(delta,s)=Lam_d(delta,s)-ed(delta,s+1).
# Lam_d(delta,s) feeds parents' ed(dp,s), dp<=delta; for dp<delta use PROVED forced F (closing holds there).
def Lam_d(delta, s, ed_self):
    tot = sp.Integer(0)
    for dp in range(0, delta+1):
        if dp < delta:
            fdp = T4.forced(dp, s)
        else:
            fdp = ed_self  # ed(delta,s), supplied
        fA = fdp.subs({d: d+1, m: m+1}, simultaneous=True)
        fB = fdp.subs({d: d+1, m: m}, simultaneous=True)
        tot += i*(2*d-1)/(2*d-dp-1) * (tau1(delta-dp)*fA - tau2(delta-dp)*fB)
    return tot

def e_desc(delta, s):
    # returns ed(delta,s), descending from top
    ed = {}
    top = delta//2
    ed[top] = sp.cancel(T2.e_top(delta))
    for ss in range(top-1, -1, -1):
        pass
    # need ed(delta,s+1) to compute ed(delta,s); but Lam_d(delta,s) needs ed(delta,s) itself. Solve self-consistently?
    return ed

# The self-reference: Lam_d(delta,s) contains ed(delta,s) at parents (dp=delta term, coeff tau(0)).
# So ed(delta,s) = Lam_d(delta,s)[with ed(delta,s)] - ed(delta,s+1): linear in ed(delta,s). Solve.
psym = sp.Symbol('ES')
for delta in [7]:
    top = delta//2
    ed = {top: sp.cancel(T2.e_top(delta))}
    for ss in range(top-1, -1, -1):
        # ed[ss] = Lam_d(delta,ss, ed_self=psym)[coeff*psym + const] - ed[ss+1]; solve ed[ss]
        Lself = Lam_d(delta, ss, psym)
        Lself = sp.expand(Lself)
        c = sp.expand(Lself).coeff(psym, 1)
        const = sp.expand(Lself - c*psym)
        # ed[ss] = c*ed[ss] + const - ed[ss+1]  => ed[ss](1-c) = const - ed[ss+1]
        ed[ss] = sp.cancel((const - ed[ss+1])/(1-c))
    diff = sp.cancel(sp.together(ed[0] - T4.arho(delta)))
    print(f"V({delta}): descending-from-top e_desc({delta},0) == a_{delta} : {diff==0}")
