
import flint, inspect
print(flint.__version__ if hasattr(flint,'__version__') else '?')
from flint import arb
try:
    print(inspect.signature(arb.__init__))
except Exception as e:
    print("sig err", e)
print(arb.__doc__)
print("arb(1.205, 0.1224) =", arb(1.205, 0.1224))
print("arb(1.205, rad=0.1224) =", arb(1.205, rad=0.1224))
print("arb('1.205', '0.1224') =", arb("1.205", "0.1224"))
