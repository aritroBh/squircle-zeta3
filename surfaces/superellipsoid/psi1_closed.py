# psi1_closed.py -- closed form of Psi1(c) = int_0^1 [sqrt(u^2+c^2(1-u)^2) - c(1-u)] log(u) du/u
# (derived with the Euler substitution N = c - t u; beta = sqrt(1+c^2)).
import mpmath as mp
def Psi1_closed(c):
    c = mp.mpf(c)
    if c == 0: return mp.mpf(-1)
    be = mp.sqrt(1 + c*c); Li2 = lambda x: mp.polylog(2, x); lg = mp.log
    Jmc = -Li2(-(be + c))
    Jpc = -Li2(be - c)
    Jmm = lg(1 + be + c)**2/2
    Jpp = lg(1 - be + c)**2/2
    Jmp = lg(2*be*(be + c))*lg((be + c - 1)/(be + c)) - Li2((be + c - 1)/(2*be)) + Li2((be + c)/(2*be))
    Jpm = lg(2*be/(be + c))*lg((be - c + 1)/(be - c)) - Li2((be - c + 1)/(2*be)) + Li2((be - c)/(2*be))
    I_minus = Jmc - Jmm - Jmp - Jpc + Jpm + Jpp
    I_plus = Jmc - Jmm - Jmp + Jpc - Jpm - Jpp
    E1 = 1 - c*c/be*lg((be + c)*(be - c + 1)/((be - c)*(be + c - 1))) + c*lg(2*c)
    return -E1 + c*c/be*I_minus - c*I_plus
if __name__ == "__main__":
    mp.mp.dps = 40
    import importlib
    P = importlib.import_module("14_Q1_pieces").Psi1
    for c in ['0.01', '0.3', '0.7071067811865475244', '0.9', '1']:
        print(c, mp.nstr(Psi1_closed(c), 30), mp.nstr(Psi1_closed(c) - P(c), 3))
