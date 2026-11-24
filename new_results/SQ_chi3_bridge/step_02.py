# Let me start by importing the engine and checking available functionality
import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C

print("Engine imported successfully")
# Check what functions/classes are available
print("Dir of cmzv module:")
print([x for x in dir(C) if not x.startswith('_')])