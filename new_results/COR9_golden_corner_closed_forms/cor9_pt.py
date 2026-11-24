"""cor9_pt.py -- QuadField patched for arbitrary D: CF-based fundamental unit and
diop_DN-based norm-l elements (works for any regulator / any class number)."""
import sys
from fractions import Fraction as Fr
import mpmath as mp
BASE = '<project>/research/'
sys.path.insert(0, BASE + '08_theory/COR7_classical_corners')
sys.path.insert(0, BASE + '08_theory/COR6_true_point_symbols')
sys.path.insert(0, BASE + '08_theory/COR9_golden_corner_closed_forms')
from cor7_pt import QuadField, QuadPoint
from cor9_lib import fund_unit, norm_l_element

class CFQuadField(QuadField):
    def __init__(self, D, verbose=False):
        self.D = D
        from cor6_nf import NF
        self.F = NF([-D, 0, 1])
        mp.mp.dps = 80
        self.rD = mp.sqrt(D)
        self.eps = fund_unit(D)
        self.gen_primes = {}
        self.mode = {}
        self._nlcache = {}
        if verbose:
            print(f"Q(sqrt{D}): eps = {self.eps}")
    def _find_norm_l(self, l, bound=200):
        if l not in self._nlcache:
            self._nlcache[l] = norm_l_element(self.D, l)[0]
        return self._nlcache[l]
