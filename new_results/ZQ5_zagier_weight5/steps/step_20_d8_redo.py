
import sys; sys.path.insert(0,'.')
import json, time
import mpmath as mp
import w5lib as W
import pipe
D=8
C = W.letters_info(D)
J, cands = pipe.load_cands(D)
K = [list(map(int, r)) for r in json.load(open(pipe.cand_dir(D)+"/kernel.json"))["K"]]
vals = pipe.phase_numerics(D, C, cands, P=1150)   # recompute with fixed cache
v2 = pipe.load_vals(D, 1150)
mp.mp.dps = 1200
print("cache roundtrip check: |cached-fresh Li5(1/2)| =",
      mp.nstr(abs(v2["li5:7"] - W.li5_series(cands[7]['x'].emb(1200), 1150)), 3))
for P in (520, 830):
    for ua in (True, False):
        R = pipe.phase_lll(D, C, cands, K, v2, P=P, use_atoms=ua)
        best = R["rels"][0]
        print(f"P={P} atoms={ua}: rels={len(R['rels'])} zeros={len(R['zeros'])} "
              f"best a0={best[1][0]} max|a|={best[0]} resid={mp.nstr(best[2],2)}")
