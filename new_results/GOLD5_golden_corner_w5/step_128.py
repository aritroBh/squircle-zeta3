from fractions import Fraction as Fr
def lyndon_words(alphabet, n):
    """all Lyndon words of length n over ordered alphabet (duval, clean)."""
    a = list(alphabet)
    words = []
    w = [-1]
    while w:
        for i in range(len(w), n):
            w.append(w[i - len(w)])
        if len(w) == n and n % len(w) == 0:
            words.append(tuple(a[i] for i in w))
        while w and w[-1] == len(a) - 1:
            w.pop()
        if w:
            w[-1] += 1
    return words
for n in range(1, 6):
    ws = lyndon_words([0, 1], n)
    print(n, [''.join(map(str, w)) for w in ws], len(ws))
