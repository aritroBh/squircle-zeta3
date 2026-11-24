
s = open('step_09b.py').read()
s = s.replace("bracket = logp - t*s2*(lv+lw) +", "bracket = logp + t*s2*(lv+lw) +")
s = s.replace("open('scratchpad/circle_bivariate.pkl','wb')","open('scratchpad/circle_bivariate.pkl','wb')")
# add consistency checks at end
s += r'''
# checks: swap symmetry b_jk = b_kj; self-dual line s2=-s1 gives a_n; area line s2=s1
import sympy as sp
pi, z3 = sp.pi, sp.zeta(3)
print("\n== swap symmetry check ==")
sw = all(sp.simplify(b[(j,k)]-b[(k,j)])==0 for j in range(6) for k in range(6) if (j,k) in b and (k,j) in b)
print("b_jk == b_kj for all:", sw)
print("\n== self-dual line s2=-s1 ==")
a2ex = 3*pi*(pi**2-8)/4
a4ex = pi*(15*pi**4 - 352*pi**2 + 1664)/192
for n,ref in [(1,0),(2,a2ex),(3,0),(4,a4ex),(5,0)]:
    val = sum((-1)**(n-j)*b[(j,n-j)] for j in range(n+1))
    print(f"a_{n} from bivariate: {sp.simplify(val)}  vs selfdual {sp.N(ref,15)} ; diff={sp.simplify(val-ref)}")
print("\n== area line s2=s1: should be area(K_p) = 4 G(1+1/p)^2/G(1+2/p) ==")
p = 2/(1+2*s1)
area = 4*sp.gamma(1+1/p)**2/sp.gamma(1+2/p)
aser = sp.series(area, s1, 0, 5).removeO()
for n in range(1,5):
    val = sum(b[(j,n-j)] for j in range(n+1))
    diff = sp.simplify(val - aser.coeff(s1,n))
    print(f"order {n}: bivariate sum vs area-series diff = {diff}")
'''
open('step_09c.py','w').write(s.replace("step_09b","step_09c"))
import subprocess
r = subprocess.run(['<project>/.venv/bin/python','step_09c.py'],capture_output=True,text=True,timeout=1500)
out = r.stdout
# print only the checks and a few b values
i = out.find("== swap")
print(out[:out.find("b_(1, 1)")])
print(out[i:])
print(r.stderr[-500:])
