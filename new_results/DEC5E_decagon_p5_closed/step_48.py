import sys, pickle, itertools
sys.path.insert(0, '../DEC5D_decagon_order5_depth')
import dec5d_lib as D
from fractions import Fraction as F
QS5 = D.QS5
bs = pickle.load(open('../DEC5D_decagon_order5_depth/dec_block_symbols.pkl','rb'))
tot = bs['tot']
# 1. symmetry check
sym_ok = True
for w in tot:
    for p in set(itertools.permutations(w)):
        if p not in tot or tot[p] != tot[w]:
            sym_ok = False; print("ASYMMETRIC:", w, p); break
    if not sym_ok: break
print("S(Wt5) fully symmetric (pure log product):", sym_ok)
# distinct multiset classes
classes = {}
for w, c in tot.items():
    classes.setdefault(tuple(sorted(w)), []).append((w, c))
print("multiset classes:", {k: len(v) for k, v in classes.items()})
# 2. solve: S(l2^a lf^b) = sum over multiset permutations of the word with a zeros, (5-a) twos,
#    each with multiplicity 1 (shuffle of X^a Y^b = binomial*word... each permutation appears once)
# Actually shuffle product X^a * Y^b as symbols: (a+b)!/(a!b!)... no: the shuffle product of
# tensors X^{*a} and Y^{*b}: X^{*a} = a! * (sum over the single word (0^a)) -> symbol of log2^a = a! (0,...,0).
# Then shuffle a!*(0^a) with b!*(2^b): gives a!b! * sum over distinct permutations. So
# S(l2^a lf^b) = a! b! * sum_{distinct perms of (0^a,2^b)} word.
# Solve for c_a: tot = sum_a c_a * a! (5-a)! * sum_perms(word). By symmetry, for each class:
#   class_value = c_a * a! (5-a)!   (each word in the class carries the same coeff)
# => c_a = tot[any word with a zeros] / (a! (5-a)!)
import math
for a in range(6):
    w = tuple([0]*a + [2]*(5-a))
    c = tot[w] / QS5(math.factorial(a)*math.factorial(5-a))
    print(f"c(l2^{a} lf^{5-a}) = {c}")
