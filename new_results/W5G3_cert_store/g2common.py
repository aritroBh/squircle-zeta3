"""W5G3 shared helpers (copy of W5G2 g2common with HERE -> W5G3_cert_store; SRC = W5G2 folder, read only)."""
import sys, os, pickle, time
sys.dont_write_bytecode = True
sys.set_int_max_str_digits(0)
from fractions import Fraction as F
ROOT = '<project>'
T = ROOT + '/research/08_theory'
SRC = T + '/W5G2_frame_exact'     # read-only: trim_*.pkl, quot_*.npz, sel.pkl, cert_summary.pkl of the first run
HERE = T + '/W5G3_cert_store'
W5C1 = T + '/W5C1_certificate_scoping'
W5G1 = T + '/W5G1_frame_galois'
A3H4 = T + '/A3h4_frame_ds'
for pth in (HERE, T + '/A3c_claude_elementary/Q', ROOT + '/paper/squircle_note/scripts', A3H4):
    if pth not in sys.path: sys.path.insert(0, pth)
from relations import words, mul, primitive
from cmzv import sh, add
E1 = {(2,): F(1), (6,): F(-1)}            # e = H(2)-H(6), value i*pi/2 (NOT a relation)
SEC = {'pp': (1, 1), 'pm': (1, -1), 'mp': (-1, 1), 'mm': (-1, -1)}
P1 = 16777213
# the 37 frame-w5 rows (indices into validated_w5.pkl) whose u=3,5 images are open (W5G1 result.md)
OPEN_ROWS = [4, 5, 7, 9, 10, 12, 14, 15, 16, 17, 52, 53, 54, 56, 58, 59, 61, 63, 64, 65, 66, 67, 68, 69, 70,
             71, 72, 73, 74, 75, 87, 91, 95, 103, 104, 107, 113]
TARGETS = [(i, u) for i in OPEN_ROWS for u in (3, 5)]      # the 74 vectors u.f_i
def act(w, u): return tuple(-1 if x == -1 else (u * x) % 8 for x in w)
def actv(d, u): return {act(w, u): c for w, c in d.items()}
def log(s, f):
    s = time.strftime('%Y-%m-%d %H:%M:%S ') + s
    print(s, flush=True)
    with open(os.path.join(HERE, f), 'a') as fh: fh.write(s + '\n')
_W5 = None
def W5():
    global _W5
    if _W5 is None: _W5 = words(5)
    return _W5
def load_framew5():
    V = pickle.load(open(W5C1 + '/validated_w5.pkl', 'rb'))
    return V
def target_vec(VW5, i, u):
    """u.f_i as dict word-tuple -> int"""
    return actv({tuple(w): c for w, c in VW5['rows'][i]}, u)
