code = r'''
from squircle import setdps
from Uq_direct import Uq_direct
from mpmath import mp, nstr, mpf
setdps(60)
T, info = Uq_direct(128, 4, 1, K=768, R=10, info=True)
print("Uq_direct T:", nstr(T, 25))
print("P:", nstr(info['P'], 25))
print("tail:", nstr(info['tail'], 25))
print("I0q:", nstr(info['I0q'], 25))
print("correct P+tail:", "0.000001035089208691401485395462")
'''
exec(code)
