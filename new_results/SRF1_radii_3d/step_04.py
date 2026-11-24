
import os
print("CWD:", os.getcwd())
d16dir = '<project>/research/08_theory/D16_pure_weight'
src = open(os.path.join(d16dir, 'd16_engine.py')).read()
print(src[:6000])
