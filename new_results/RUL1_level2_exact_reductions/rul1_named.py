"""rul1_named.py -- named basis constants as exact word-combos (G-values)."""
from fractions import Fraction as Fr
import rul1

def _sh(u, v):
    return rul1.shuffle(u, v)

def _scale(A, c):
    return {k: v*c for k, v in A.items()}

def _shpow(w, n):
    """shuffle power w^{ш n}"""
    out = {w: Fr(1)}
    for _ in range(n-1):
        out = rul1.sh_mul(out, {w: Fr(1)})
    return out

def _li_half(n):
    """Li_n(1/2) = sum_{c in {0,-1}^{n-1}} (-1)^{#(-1)} G(c, -1; 1)."""
    import itertools
    out = {}
    for c in itertools.product((0, -1), repeat=n-1):
        w = c + (-1,)
        sgn = Fr((-1)**sum(1 for a in c if a == -1))
        out[w] = out.get(w, Fr(0)) + sgn
    return out

_m1 = (-1,)
_o1 = (0, 1)
_oo1 = (0, 0, 1)
_ooo1 = (0, 0, 0, 1)
_oooo1 = (0, 0, 0, 0, 1)

NAMED = {
 1: [('log2', {_m1: Fr(1)})],
 2: [('pi2', _scale(_sh(_o1, _o1), Fr(-6, 2))),   # pi^2 = -6 G(0,1);  G(0,1)шG(0,1)=2G(0,1)... careful: pi2 = -6*G(0,1) as single word
     ('log2^2', _shpow(_m1, 2))],
 3: [('z3', {(0, 0, 1): Fr(-1)}),
     ('pi2*log2', _scale(rul1.sh_mul({_o1: Fr(1)}, {_m1: Fr(1)}), Fr(-6))),
     ('log2^3', _shpow(_m1, 3))],
 4: [('Li4h', _li_half(4)),
     ('z3*log2', _scale(rul1.sh_mul({_oo1: Fr(1)}, {_m1: Fr(1)}), Fr(-1))),
     ('pi4', _scale(rul1.sh_mul({_o1: Fr(1)}, {_o1: Fr(1)}), Fr(36))),
     ('pi2*log2^2', _scale(rul1.sh_mul({_o1: Fr(1)}, _shpow(_m1, 2)), Fr(-6))),
     ('log2^4', _shpow(_m1, 4))],
 5: [('Li5h', _li_half(5)),
     ('Li4h*log2', rul1.sh_mul(_li_half(4), {_m1: Fr(1)})),
     ('z5', {_oooo1: Fr(-1)}),
     ('z3*pi2', _scale(rul1.sh_mul({_oo1: Fr(1)}, {_o1: Fr(1)}), Fr(6))),
     ('z3*log2^2', _scale(rul1.sh_mul({_oo1: Fr(1)}, _shpow(_m1, 2)), Fr(-1))),
     ('pi4*log2', _scale(rul1.sh_mul(rul1.sh_mul({_o1: Fr(1)}, {_o1: Fr(1)}), {_m1: Fr(1)}), Fr(36))),
     ('pi2*log2^3', _scale(rul1.sh_mul({_o1: Fr(1)}, _shpow(_m1, 3)), Fr(-6))),
     ('log2^5', _shpow(_m1, 5))],
}
# fix pi2: must be -6*G(0,1) (single word), not a shuffle square
NAMED[2] = [('pi2', {_o1: Fr(-6)}), ('log2^2', _shpow(_m1, 2))]

NAMED_BY_WEIGHT = NAMED

WEIGHT_OF = {}
for w, lst in NAMED.items():
    for nm, _ in lst:
        WEIGHT_OF[nm] = w
WEIGHT_OF['1'] = 0
NAMED['1'] = [('1', {(): Fr(1)})]
NAMED_BY_WEIGHT[0] = NAMED['1']

# E5e 20-constant basis order (weights 0..5)
E5E_ORDER = ['1', 'log2', 'pi2', 'log2^2', 'z3', 'pi2*log2', 'log2^3',
             'pi4', 'z3*log2', 'pi2*log2^2', 'log2^4', 'Li4h',
             'z5', 'z2z3', 'pi4*log2', 'z3*log2^2', 'pi2*log2^3', 'log2^5', 'Li5h', 'Li4h*log2']
# z2z3 = zeta(2)zeta(3) = pi2*z3/6 = (z3*pi2)/6
