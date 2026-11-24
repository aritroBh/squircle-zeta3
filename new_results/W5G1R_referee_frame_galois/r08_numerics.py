"""r08: independent double-precision numerics (rcommon.H; ~1e-13 relative) for the residual X and the 74 u.f,
with a control (the u=3 image of one frame row with the sign of the u=5 part flipped is not expected to vanish)."""
from rcommon import *
import pickle
sys.path.insert(0, ROOT + '/paper/squircle_note/scripts')
L = 'r08.log'
VW5 = pickle.load(open(T + '/W5C1_certificate_scoping/validated_w5.pkl', 'rb'))
FR = [{tuple(w): F(c) for w, c in row} for row in VW5['rows']]
R = pickle.load(open(HERE + '/r05_X.pkl', 'rb'))
def scaled(d):
    num = sum(float(c) * H(w) for w, c in d.items())
    den = sum(abs(float(c) * H(w)) for w, c in d.items())
    return abs(num) / den
worst = max(scaled(actv(FR[i], u)) for i, u in R['nonlit'])
log(f'74 non-literal u.f: max |u.f.V| / sum|c H| = {worst:.2e}', L)
for c in ('rat', 'sq'):
    x = R['X'][c]; mx = max(abs(v) for v in x.values())
    xs = {w: v / mx for w, v in x.items()}
    log(f'X_{c}: |X.V| / sum|c H| = {scaled(xs):.2e}', L)
# controls: a frame row image under a non-Galois letter map (j -> j+1 on nonzero letters, 0 fixed) and an e-shifted image
i, u = R['nonlit'][0]
bad = {tuple(x if x in (-1, 0) else (x % 7) + 1 for x in w): c for w, c in actv(FR[i], u).items()}
log(f'control (letters 1..7 cycled): scaled value {scaled(bad):.2e}', L)
