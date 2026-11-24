import pickle, sys
sys.path.insert(0,'.')
import engine14 as E
S = pickle.load(open('S14.pkl','rb'))
from collections import Counter
# classify monomials by the "atom" content (non-letter, non-z symbols)
def classify(mon):
    d = dict(mon)
    atoms = [s for s in d if s not in ('l2','la','lu','z2','z3')]
    letters = {s:d[s] for s in d if s in ('l2','la','lu')}
    zs = {s:d[s] for s in d if s in ('z2','z3')}
    return (tuple(sorted(atoms)), tuple(sorted(letters.items())), tuple(sorted(zs.items())))
atomcount = Counter()
pureletter = 0
for mon in S:
    a,l,z = classify(mon)
    if not a: pureletter += 1
    else:
        for s in a: atomcount[s]+=1
print("total terms:", len(S), " pure letter/z terms:", pureletter)
print("atoms appearing:", dict(atomcount))
# Group terms by atom signature (ignoring letters/z)
byatom = Counter()
for mon in S:
    d=dict(mon); atoms=tuple(sorted(s for s in d if s not in ('l2','la','lu','z2','z3')))
    byatom[atoms]+=1
print("\natom signatures (top):")
for sig,c in byatom.most_common(40):
    print("  ", sig, c)
