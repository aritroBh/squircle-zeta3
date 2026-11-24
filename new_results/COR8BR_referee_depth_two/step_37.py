# Take actual pairs from COR8B's database; compare my delta vs their engine vs their STATED closed forms
import sys, pickle
sys.path.insert(0,'<project>/research/08_theory/COR8B_depth_two_two_terms')
sys.path.insert(0,'.')
import cor8b_sym as CB
import cor8br_sym as E
from fractions import Fraction as Fr
import sympy as sp

data = pickle.load(open('<project>/research/08_theory/COR8B_depth_two_two_terms/pairs_delta_all_v2.pkl','rb'))
print("pairs in v2 db:", len(data))
print("first entry:", data[0])
