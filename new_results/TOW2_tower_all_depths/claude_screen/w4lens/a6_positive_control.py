# Positive control: the +-i-simple-pole / real-residue detectors in a5 are not vacuous.
import sys, os, pickle, copy
from fractions import Fraction as Fr
sys.dont_write_bytecode = True
exec(open('a5_image_and_blindness.py').read().split("print(\"(A)")[0])   # reuse definitions only
T5=copy.deepcopy(Tw[5])
for r in (2,3,4,5,6):
    for comp in range(4):
        T2=copy.deepcopy(Tw[5])
        base=list(T2[(4,0)].get(r,(0,0,0,0)))
        base[comp]=Fr(base[comp])+Fr(1,7)
        T2[(4,0)][r]=tuple(base)
        P=kernel_from_tower(T2,(3,0))
        rd,simple,const=realdict(P)
        want={rr:c for rr,c in Kw[5][(3,0)]['real'].items() if any(x!=0 for x in c)}
        has_i=any(not e.is_zero() for e in simple[1]) or any(not e.is_zero() for e in simple[2]) or (not const.is_zero())
        print(f"perturb (4,0) tower r={r} comp={comp}: kernel(3,0) real residues changed={rd!=want}, +-i simple/const present={has_i}")
