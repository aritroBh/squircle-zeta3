
import sys; sys.path.insert(0,'.')
import w5lib as W
import pipe
C = W.letters_info(13)
for B in (10, 14):
    cands = W.enum_exceptional(C, B=B)
    print(f"B={B}: {len(cands)} exceptional S-units")
    if len(cands) > 4:
        for c in cands:
            print(f"   x={pipe.fmt_qd(c['x']):>22s} vx={c['vx']} 1-x={c['v1m']} h={c['h']}")
