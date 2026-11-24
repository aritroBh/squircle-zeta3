
import pickle, sys
sys.dont_write_bytecode = True
import sympy as sp
from fractions import Fraction as Fr

# heights of w7 certificates
for m in range(3):
    cert = pickle.load(open(f'cert_D7_{m}.pkl','rb'))
    hs = [sp.Rational(c) for c in cert.values()]
    bits = sorted(max(abs(h.p), h.q).bit_length() for h in hs)
    print(f"D~_{{7,{m}}}: {len(hs)} entries, height bits: median {bits[len(bits)//2]}, max {bits[-1]}")
