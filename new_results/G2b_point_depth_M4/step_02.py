"""step_02.py -- G2b Task 2: TRUE point-level (3,2) pairings at c0 = 1/2 (hexagon) and c0 = 1/3.

Same fixed machinery as step_01 (expand_sym before pairing). F = Q, S = {2,3}, basis {2,3},
X = log 2, Y = log 3. Also the NAIVE (function-level symbol specialised at the point) pairing,
fixed the same way, for contrast. Zero/nonzero decides classicality at the symbol level.
"""
import sys, os, json, time
from fractions import Fraction as Fr
BASE = "<project>/research/08_theory"
G2B = BASE + "/G2b_point_depth_M4/"
sys.path.insert(0, BASE + "/COR6_true_point_symbols")
sys.path.insert(0, BASE + "/G2_M4alpha")
os.chdir(BASE + "/G2_M4alpha")
import importlib.util
_spec = importlib.util.spec_from_file_location("g2b_step01", G2B + "step_01.py")
_s01 = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(_s01)
ft_expand_sym = _s01.ft_expand_sym; obstruction_fixed = _s01.obstruction_fixed
shuffle_control = _s01.shuffle_control; li5_control = _s01.li5_control
li5_sunit_control = _s01.li5_sunit_control; log = _s01.log
from step_12 import br, S5_func
from step_13 import W5Point

t0 = time.time()

def naive_obstruction_fixed(pt):
    """Specialise the FUNCTION-level weight-5 symbol (step_12) at the point's letter vectors,
    then expand to basis letters and pair (the same expand fix)."""
    lv = pt.letter_vectors()
    M = [lv["c"], lv["1-c"], lv["1+c"]]
    S5 = S5_func()                            # {word over {1,2,3}: {e: Fr}}
    S = {}
    for word, d in S5.items():
        val = pt.F.zero()
        for e, cf in d.items():
            val = pt.F.add(val, pt.F.scale(cf, pt.F.pow(pt.c0, e)))
        S[tuple(M[ix - 1] for ix in word)] = val
    SE = ft_expand_sym(S, len(pt.basis))
    return obstruction_fixed(SE), len(SE)

def run_point(label, minpoly, c0, basis=None, zargs=()):
    pt = W5Point(label, minpoly, c0, basis=basis)
    log(f"{label}: field deg {pt.F.d}, S={pt.S0}, dimV={pt.dimV}, basis={pt.basis} ({time.time()-t0:.0f}s)")
    S = pt.omega5()
    log(f"  true point weight-5 symbol: {len(S)} tensor words")
    SE = ft_expand_sym(S, pt.dimV)
    log(f"  expanded over basis letters: {len(SE)} words")
    obs = obstruction_fixed(SE)
    log(f"  TRUE (3,2) pairing [obs_Q, obs_T]: {[[str(x) for x in oo] for oo in obs]}")
    nz = any(any(x for x in oo) for oo in obs)
    log(f"  -> pairing is {'NONZERO (depth-two at the symbol level)' if nz else 'ZERO (classical at the symbol level)'}")
    log(f"  letter vectors: { {k: [str(x) for x in v] for k, v in pt.letter_vectors().items()} }")
    nobs, nse = naive_obstruction_fixed(pt)
    log(f"  NAIVE (specialised function, fixed) pairing ({nse} words): {[[str(x) for x in oo] for oo in nobs]}")
    if zargs:
        log(f"  Li5-at-S-units control: {li5_sunit_control(pt, zargs)}")
    return pt, obs, nobs

if __name__ == '__main__':
    log(f"step_02 start {time.strftime('%H:%M:%S')}")
    # c0 = 1/2 (hexagon corner, alpha = 2pi/3), F = Q, basis {2, 3}
    P6, obs6, nobs6 = run_point("c=1/2 (hexagon)", [-1, 1], (Fr(1, 2),),
                                zargs=[(Fr(1, 2),), (Fr(1, 3),), (Fr(2, 3),), (Fr(-1, 2),), (Fr(-2, 1),)])
    # c0 = 1/3, F = Q, basis {2, 3}
    P3, obs3, nobs3 = run_point("c=1/3", [-1, 1], (Fr(1, 3),),
                                zargs=[(Fr(1, 2),), (Fr(1, 3),), (Fr(2, 3),), (Fr(-1, 2),), (Fr(-2, 1),)])
    log(f"  shuffle control (200 random products): {shuffle_control()}")
    log(f"  Li5 control (200 random Li5-form symbols): {li5_control()}")
    json.dump({'hexagon': {'obs': [[str(x) for x in oo] for oo in obs6],
                           'naive': [[str(x) for x in oo] for oo in nobs6]},
               'third': {'obs': [[str(x) for x in oo] for oo in obs3],
                         'naive': [[str(x) for x in oo] for oo in nobs3]}},
              open(G2B + 'pair_half_third.json', 'w'), indent=1)
    log(f"step_02 done ({time.time()-t0:.0f}s)")
