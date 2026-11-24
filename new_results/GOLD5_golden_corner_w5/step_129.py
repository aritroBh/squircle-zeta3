import itertools
def is_lyndon(w):
    return all(w < w[i:] + w[:i] for i in range(1, len(w)))
def lyndon_words(n):
    return [w for w in itertools.product((0, 1), repeat=n) if is_lyndon(w)]
for n in range(1, 6):
    ws = lyndon_words(n)
    print(n, [''.join(map(str, w)) for w in ws], len(ws))
