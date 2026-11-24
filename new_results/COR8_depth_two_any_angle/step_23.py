
step01 = '''"""step_01.py -- build and validate the eps-symbol engine (TRUE symbols of Li31/Li22)."""
import sys, json
sys.path.insert(0, "<project>/research/08_theory/COR8_depth_two_any_angle")
from fractions import Fraction as Fr
from cor8_sym import *

ok = True
def check(name, cond):
    global ok
    print(("PASS " if cond else "FAIL ") + name)
    ok = ok and cond

Z5 = (Fr(0),) * 5
def v5(j, val=1):
    v = list(Z5); v[j] = Fr(val); return tuple(v)
# formal letters: 0=X, 1=Y, 2=P=[1-x], 3=Q=[1-y], 4=R=[1-xy]
X, Y, P, Q, R = (v5(j) for j in range(5))

# ---- 1. S(Li31(x,1)) = P P X X : substitute Y->0, R->P; Q-terms must cancel
S = li31_sym(X, Y, P, Q, R)
zero = v5(0, 0)
S1 = {}
for w, cf in S.items():
    w2 = tuple(P if l == R else l for l in w)
    if Y in w2:
        continue
    S1[w2] = padd(S1.get(w2, {}), cf)
S1 = {w: cf for w, cf in S1.items() if cf}
check("Li31(x,1): Q-terms cancel", all(Q not in w for w in S1))
check("Li31(x,1) == PPXX", S1 == {(P, P, X, X): ONE})

# ---- 2. S(Li22(x,1)) = P X P X
S = li22_sym(X, Y, P, Q, R)
S1 = {}
for w, cf in S.items():
    w2 = tuple(P if l == R else l for l in w)
    if Y in w2 or Q in w2:
        continue
    S1[w2] = padd(S1.get(w2, {}), cf)
S1 = {w: cf for w, cf in S1.items() if cf}
# Q-terms: P Q Y X, -R Q Y X, P Q X Y all contain Y -> dropped. Verify nothing with Q remains silently
check("Li22(x,1) == PXPX", S1 == {(P, X, P, X): ONE})

# ---- 3. S(Li31(x,-1)) = [1+x][1-x][x][x] : Y->0, Q->0, R->M
M = v5(1)  # reuse index 1 as M=[1+x] in this check
S = li31_sym(X, Y, P, Q, R)
S1 = {}
for w, cf in S.items():
    w2 = tuple(M if l == R else l for l in w)
    if Y in w2 or Q in w2:
        continue
    S1[w2] = padd(S1.get(w2, {}), cf)
S1 = {w: cf for w, cf in S1.items() if cf}
check("Li31(x,-1) == M P X X", S1 == {(M, P, X, X): ONE})

# ---- 4. S(Li22(x,-1)) = M X P X
S = li22_sym(X, Y, P, Q, R)
S1 = {}
for w, cf in S.items():
    w2 = tuple(M if l == R else l for l in w)
    if Y in w2 or Q in w2:
        continue
    S1[w2] = padd(S1.get(w2, {}), cf)
S1 = {w: cf for w, cf in S1.items() if cf}
check("Li22(x,-1) == M X P X", S1 == {(M, X, P, X): ONE})

# ---- 5. stuffle: S11(x,y)+S11(y,x) = PQ + QP + R(X+Y)
lhs = sadd(li11_sym(X, Y, P, Q, R), li11_sym(Y, X, Q, P, R))
rhs = {(P, Q): dict(ONE), (Q, P): dict(ONE), (R, X): dict(ONE), (R, Y): dict(ONE)}
check("stuffle Li11", lhs == rhs)

# ---- 6. delta22(Li4) = 0 for three S-unit args (letters over {c,m,p}=0,1,2)
Z3 = (Fr(0),) * 3
def v3(j, val=1):
    v = list(Z3); v[j] = Fr(val); return tuple(v)
cv, mv, pv = v3(0), v3(1), v3(2)
for name, rv, omrv in [("g=c", cv, mv), ("g=1-c", mv, cv), ("g=(1-c)/(1+c)", ladd(mv, lneg(pv)), ladd(cv, lneg(pv)) )]:
    d = delta22(li_sym(4, rv, omrv), 3)
    check(f"delta22 Li4({name}) == 0", not d)

# ---- 7. delta22 of a shuffle product = 0
sh = sshuffle(li_sym(2, cv, mv), li_sym(2, mv, cv))
check("delta22 Li2(c) sh Li2(1-c) == 0", not delta22(sh, 3))
sh2 = sshuffle(li_sym(3, cv, mv), log_sym(pv))
check("delta22 Li3(c)*log(1+c) == 0", not delta22(sh2, 3))

# ---- 8. delta22(F) == COR4 value
FS = json.load(open("<project>/research/08_theory/COR4_M3alpha_classical/F_symbol.json"))
def parse_letter(s):
    # COR4 basis [2,c,1-c,1+c]; 4 chars
    digs = []
    i = 0
    while i < len(s):
        if s[i] == '-':
            digs.append(int(s[i:i+2])); i += 2
        else:
            digs.append(int(s[i])); i += 1
    assert len(digs) == 4, s
    assert digs[0] == 0, "letter 2 present!"
    return (Fr(digs[1]), Fr(digs[2]), Fr(digs[3]))
F = {}
for key, val in FS.items():
    e_str, word_str = key.split(" ")
    e = int(e_str[2:])
    word = tuple(parse_letter(b) for b in word_str.split("|"))
    F[word] = padd(F.get(word, {}), {e: Fr(val)})
dF = delta22(F, 3)
print("delta22(F) components:")
for k, poly in sorted(dF.items()):
    print("  ", k, {e: str(v) for e, v in sorted(poly.items())})
target = {((0, 1), (1, 2)): {0: Fr(1), -1: Fr(-1)},
          ((0, 2), (1, 2)): {0: Fr(-1), -1: Fr(-1)}}
check("delta22(F) == (1-1/c) e1^e3 + (-1-1/c) e2^e3", dF == target)

# ---- 9. delta22(Li22(c,-1)) nonzero (DGR-style depth-2)
Z3w = (Fr(0),) * 3
Xc, Pc, Mc = v3(0), v3(1), v3(2)
Y0 = v3(0, 0)
d = delta22(li22_sym(Xc, Y0, Pc, Y0, Mc), 3)
print("delta22 Li22(c,-1):", {k: {e: str(v) for e, v in sorted(p.items())} for k, p in d.items()})
check("delta22 Li22(c,-1) != 0", bool(d))

print("ALL PASS" if ok else "SOME FAILED")
'''
open('step_01.py','w').write(step01)
import subprocess
r = subprocess.run(['<project>/.venv/bin/python','step_01.py'],capture_output=True,text=True)
print(r.stdout); print(r.stderr[-2000:])
open('step_01.out','w').write(r.stdout)
