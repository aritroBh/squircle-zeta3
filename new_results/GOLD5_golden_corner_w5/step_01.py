"""step_01.py -- GOLD5 Task 2: TRUE point-level (3,2) pairings at the golden corners,
with the G2b square (84,-38 sqrt2) and hexagon (41,-33/2) CONTROLS first.

Machinery: G2b's fixed pipeline (gold5_lib: ft_expand_sym + obstruction_fixed from g2b_step01;
W5Point from G2_M4alpha/step_13). Golden field Q(sqrt5), S = {2} (2 inert), basis {2, phi},
X = log 2, Y = log phi. Letters:
  c0 = 1/phi:   c = phi^-1 (0,-1), 1-c = phi^-2 (0,-2), 1+c = phi (0,1)  (+ letter 2 from gaps)
  c0 = sqrt5-2: c = phi^-3 (0,-3), 1-c = 2 phi^-2 (1,-2), 1+c = 2 phi^-1 (1,-1)
Controls: shuffle products and Li5 symbols pair to 0; Li5 at golden S-units pair to 0.
"""
import sys, os, json, time
from fractions import Fraction as Fr
sys.path.insert(0, '<project>/research/08_theory/GOLD5_golden_corner_w5')
from gold5_lib import (GOLD, s01, golden_point, log, nf_str, W5Point)
t0 = time.time()

if __name__ == '__main__':
    log(f"step_01 start {time.strftime('%H:%M:%S')}")
    # ---- controls: square and hexagon (G2b values)
    P4 = W5Point("c=1/sqrt2 (square)", [-2, 0, 1], (Fr(0), Fr(1, 2)),
                 basis=[(Fr(2), Fr(0)), (Fr(1), Fr(1))])
    obs4 = s01.obstruction_fixed(s01.ft_expand_sym(P4.omega5(), 2))
    ok4 = obs4[0] == (Fr(84), Fr(0)) and obs4[1] == (Fr(0), Fr(-38))
    log(f"CONTROL square: {[[str(x) for x in oo] for oo in obs4]} target [(84,0),(0,-38)] "
        f"{'PASS' if ok4 else 'FAIL'} ({time.time()-t0:.0f}s)")
    P6 = W5Point("c=1/2 (hexagon)", [-1, 1], (Fr(1, 2),))
    obs6 = s01.obstruction_fixed(s01.ft_expand_sym(P6.omega5(), 2))
    ok6 = obs6[0] == (Fr(41),) and obs6[1] == (Fr(-33, 2),)
    log(f"CONTROL hexagon: {[[str(x) for x in oo] for oo in obs6]} target [(41,),(-33/2,)] "
        f"{'PASS' if ok6 else 'FAIL'} ({time.time()-t0:.0f}s)")
    assert ok4 and ok6

    out = {'square_control': [[str(x) for x in oo] for oo in obs4],
           'hexagon_control': [[str(x) for x in oo] for oo in obs6]}
    # ---- golden corners
    for tag in ('gol1', 'gol3'):
        P = golden_point(tag)
        lv = P.letter_vectors()
        log(f"{P.label}: letters {{k: [exponents in (2, phi)]}} = "
            f"{ {k: [str(x) for x in v] for k, v in lv.items()} } ({time.time()-t0:.0f}s)")
        S = P.omega5()
        log(f"  true point weight-5 symbol: {len(S)} tensor words")
        SE = s01.ft_expand_sym(S, 2)
        log(f"  expanded over basis letters: {len(SE)} words")
        obs = s01.obstruction_fixed(SE)
        log(f"  TRUE (3,2) pairing [obs_Q, obs_T] = [{nf_str(obs[0])}, {nf_str(obs[1])}]"
            f"  = {[[str(x) for x in oo] for oo in obs]}")
        nz = any(any(x for x in oo) for oo in obs)
        log(f"  -> pairing is {'NONZERO (depth two needed at the symbol level)' if nz else 'ZERO'}")
        # controls at this point: shuffles and Li5-at-S-units
        zargs = [(Fr(-1, 2), Fr(1, 2)), (Fr(3, 2), Fr(-1, 2)), (Fr(1, 2), Fr(-1, 2)),
                 (Fr(1, 2), Fr(0)), (Fr(-2), Fr(1)), (Fr(9), Fr(-4)), (Fr(-1, 2), Fr(-1, 2))]
        # 1/phi, 1/phi^2, -1/phi, 1/2, sqrt5-2=phi^-3, 9-4sqrt5=phi^-6, -phi: exceptional {2}-units
        okL = s01.li5_sunit_control(P, zargs)
        log(f"  Li5-at-golden-S-units control: {okL}")
        out[tag] = {'letters': {k: [str(x) for x in v] for k, v in lv.items()},
                    'obs': [[str(x) for x in oo] for oo in obs],
                    'n_words': len(S), 'n_expanded': len(SE), 'li5_control': bool(okL)}
    log(f"  shuffle control (200): {s01.shuffle_control()}")
    log(f"  Li5 control (200): {s01.li5_control()}")
    json.dump(out, open(GOLD + 'pair_golden.json', 'w'), indent=1)
    log(f"step_01 done ({time.time()-t0:.0f}s)")
