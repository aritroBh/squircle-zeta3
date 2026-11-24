# Quadrature of  int_b^a sum_t coef_t/(v-pole_t) log^i t log^j(1-t) dv  (same scheme as s2_eval.py).
# terms: list of (i, j, pole_label, mpf_coefficient); pole_label in {"0","a","-a","b","-b"}.
# Call only after setting mp.mp.dps (AFTER all imports).
import mpmath as mp
def integrate_terms(terms):
    sq2 = mp.sqrt(2); a = 1 + sq2; b = sq2 - 1
    def make(side):
        def f(s):
            if side == 1:
                v = b + s; t = s*(v + a)/(4*v); omt = (2 - s)*(v + b)/(4*v)
                l1 = mp.log(t); l2 = mp.log1p(-t)
                dv = {"0": v, "a": s - 2, "-a": v + a, "b": s, "-b": v + b}
            else:
                v = a - s; omt = s*(v + b)/(4*v); t = (2 - s)*(v + a)/(4*v)
                l2 = mp.log(omt); l1 = mp.log1p(-omt)
                dv = {"0": v, "a": -s, "-a": v + a, "b": 2 - s, "-b": v + b}
            tot = mp.mpf(0)
            for i, j, lab, c in terms:
                tot += c/dv[lab]*l1**i*l2**j
            return tot
        return f
    I1, e1 = mp.quad(make(1), [0, 1 - b], error=True, maxdegree=12)
    I2, e2 = mp.quad(make(2), [0, a - 1], error=True, maxdegree=12)
    return I1 + I2, e1 + e2
