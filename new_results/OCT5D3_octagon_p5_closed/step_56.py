import os, sys, pickle
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
import mpmath as mp
mp.mp.dps = 1650
A = {k: mp.mpf(v) for k, v in pickle.load(open("atoms_1650.pkl","rb")).items()}
chk = A['Li4(b^2)'] - 8*(A['Li4(b)']+A['Li4(-b)'])
print("Li4(b^2)-8(Li4(b)+Li4(-b)):", mp.nstr(abs(chk), 3))
chk2 = A['Li5(b^4)'] - 256*(A['Li5(b)']+A['Li5(-b)']) - 16*A['Li5(-b^2)']
print("Li5(b^4)-256(Li5(b)+Li5(-b))-16Li5(-b^2):", mp.nstr(abs(chk2), 3))
chk3 = A['Li5(-b/2)'] - mp.polylog(5, -(mp.sqrt(2)-1)/2)
print("Li5(-b/2) self-consistency:", mp.nstr(abs(chk3), 3))
