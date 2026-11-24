"""Compare the a_j stored in store/*.core.jsonl.xz with W5G2_frame_exact/cert_summary.pkl (first run), all 74 targets, exact values."""
import lzma, json, pickle, sys
from fractions import Fraction as F
sys.set_int_max_str_digits(0)
old = pickle.load(open('<project>/research/08_theory/W5G2_frame_exact/cert_summary.pkl', 'rb'))['summary']
same = 0; bad = []
for rec in old:
    i, u = rec['i'], rec['u']
    with lzma.open(f'store/u{u}.f{i}.core.jsonl.xz', 'rt') as f:
        fr = next(o for o in map(json.loads, f) if o['section'] == 'frame')['terms']
    new = {((t[0], t[1], 1)): F(*map(int, t[-1].split('/'))) for t in fr}
    new = {(('framew5' if k[0] == 'framew5' else 'frame234'), k[1], 1): v for k, v in new.items()}
    o = {k: F(v) for k, v in rec['frame_terms']}
    if new == o: same += 1
    else: bad.append((i, u, len(new), len(o)))
print(f'a_j identical for {same}/{len(old)} targets (exact Fraction equality of the coefficient maps, keys = frame row identity); different: {bad}')
