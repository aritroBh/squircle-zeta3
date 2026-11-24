fixed = '''def dynkin_left(w):
    """delta(l1 l2 ... ln) = [[...[l1,l2],...],ln] as dict word -> Fr."""
    L = {}
    for w2, c in [(w[:2], Fr(1)), (w[1::-1], Fr(-1))]:
        L[w2] = L.get(w2, 0) + c
    L = {x: c for x, c in L.items() if c}
    for k in range(2, len(w)):
        out = {}
        for u, c in L.items():
            for w2, c2 in lieb_words(u, w[k:k+1]).items():
                out[w2] = out.get(w2, 0) + c * c2
        L = {x: c for x, c in out.items() if c}
    return L
'''
import re
src = open("cor5_sym.py").read()
# replace the old dynkin_left definition
start = src.find("def dynkin_left")
end = src.find("def proj_lie")
src = src[:start] + fixed + "\n" + src[end:]
open("cor5_sym.py", "w").write(src)
print("patched")
