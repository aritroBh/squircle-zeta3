import pickle, sys
sys.path.insert(0, ".")
rewrite = pickle.load(open("rewrite_exact.pkl","rb"))
allwords = set()
for r in rewrite.values(): allwords.update(r.keys())
print("words starting with slot 1:", sum(1 for w in allwords if w and w[0] == 1))
from dod5d_grewrite import L0, L1, P0, P1, P2, P3
# convert QS3-letter words to slot words
SLOT = {L0: 0, L1: 1, P0: 2, P1: 3, P2: 4, P3: 5}
slotwords = set()
for w in allwords:
    slotwords.add(tuple(SLOT[x] for x in w))
MREV = (1, 0, 3, 2, 5, 4)
closure = set()
for w in slotwords:
    for i in range(len(w)+1):
        closure.add(w[i:])
        closure.add(tuple(MREV[x] for x in reversed(w[:i])))
print("needed:", len(slotwords), "closure:", len(closure))
extra1 = sum(1 for w in closure if w and w[0] == 1)
print("closure words starting with 1 (from mrev):", extra1)
pickle.dump(slotwords, open("needed_words.pkl","wb"))
