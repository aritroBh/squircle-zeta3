# Let's see what atom returns for different arguments:
X = C.atom(('li', 3, C.B))
print("X (atom(('li', 3, C.B))) =")
print(X)
print("\nY = C.atom(('li', 3, -C.B))")
Y = C.atom(('li', 3, -C.B))
print(Y)

print("\nLv = C.atom(('log', C.A))")
Lv = C.atom(('log', C.A))
print(Lv)

print("\nP2 = C.atom(('pi2',))")
P2 = C.atom(('pi2',))
print(P2)