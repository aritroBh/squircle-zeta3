import sys; sys.path.insert(0,'.')
import point_engine as PE
import kimilib, pickle
from fractions import Fraction as F
forms = pickle.load(open('forms10.pkl','rb'))
a=(1,1,1,1,1)
# W=3: forms[3] has 1 coord (z3). 
ref3 = kimilib.full_coords_at(a, 3, forms)[3]
got3 = PE.reduced_at(a, 3)[3]
print("W=3 z3: ref", ref3, " got", got3)
# m0 alone
print("m0[3]:", kimilib.m0_coords(a,3)[3])
# regular word dict at a, weight 3
vw = PE.volreg_words_point(a,3)[3]
print("vw[3]:", vw)
# from engine5 symbolic for comparison
import engine5
eng = engine5.build_engine(4)
print("engine5 Areg[3]:", {w: kimilib.peval5(pp,a) for w,pp in eng['Areg'][3].items()})
print("engine5 Gg[3]:", {w: kimilib.peval5(pp,a) for w,pp in eng['Gg'][3].items()})
print("engine5 Inn[3]:", {w: kimilib.peval5(pp,a) for w,pp in eng['Inn'][3].items()})
