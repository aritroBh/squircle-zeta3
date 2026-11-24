import time
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
import sys; sys.path.insert(0, WD)
from cor9_lib import class_number, factorint
t0 = time.time()
sqfree = [D for D in range(2, 201) if all(e == 1 for e in factorint(D).values())]
h = {D: class_number(D) for D in sqfree}
h1 = [D for D in sqfree if h[D] == 1]
print(f"{len(sqfree)} squarefree D <= 200 ({time.time()-t0:.0f}s)")
print(f"class number one ({len(h1)} fields):")
print(h1)
print("\nCOR7's list of 15:", [2,3,5,6,7,11,13,17,19,21,23,29,33,37,41])
print("h=1 fields <= 50 COR7 missed:", [D for D in h1 if D <= 50 and D not in [2,3,5,6,7,11,13,17,19,21,23,29,33,37,41]])
open(WD+'progress.log','a').write(f"{time.strftime('%H:%M:%S')} step_55: class numbers D<=200 computed; {len(h1)} h=1 fields.\n")
