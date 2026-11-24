import os, pickle, sys
base = "<project>"
sys.path.insert(0, os.path.join(base,"research/08_theory/OCT8_symbol_test"))
import oct8
print([n for n in dir(oct8) if not n.startswith('_')][:50])
src = open(os.path.join(base,"research/08_theory/OCT8_symbol_test/oct8.py")).read()
print(src[:4000])
