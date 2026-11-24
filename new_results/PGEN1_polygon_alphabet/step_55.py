import pickle
from fractions import Fraction as Fr
CR12 = pickle.load(open("CR_regular.pkl","rb"))["CR12"]
for z in sorted(CR12, key=repr):
    print(repr(z)[:90], [type(t).__name__ for t in z])
