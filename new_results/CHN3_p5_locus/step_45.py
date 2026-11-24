import time
print(open('engine12.log').read())
import os
print("forms12.pkl exists:", os.path.exists('forms12.pkl'))
# peek at engine.py for P4 reduced_coeffs signature
src = open('engine.py').read()
print(src[:2000])
