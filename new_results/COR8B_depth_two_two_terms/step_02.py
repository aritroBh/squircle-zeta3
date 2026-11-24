
script = r'''
# step_01.py -- verify the G_phys word <-> convB MPL identification numerically
import mpmath as mp
mp.mp.dps = 50

def li_convB(a, b, x, y):
    # Sum_{n>m>=1} x^n y^m / (n^a m^b), converges for |x|<1, |x*y|<1
    # sum over m, inner tail sum over n>m done analytically:
    # S = sum_m y^m/m^b * sum_{n>m} x^n/n^a ; tail T_m = Li_a(x) - sum_{n<=m} x^n/n^a
    x = mp.mpf(x); y = mp.mpf(y)
    total = mp.mpf(0)
    m = 0
    partial = mp.mpf(0)  # sum_{n<=m} x^n/n^a
    xp = mp.mpf(1)       # x^m
    while True:
        m += 1
        xp *= x
        partial += xp / mp.power(m, a)
        tail = mp.polylog(a, x) - partial
        term = mp.power(y, m) / mp.power(m, b) * tail
        total += term
        if m > 30 and abs(term) < mp.mpf(10)**(-mp.mp.dps+10) and abs(tail) < 1:
            break
    return total

def li1(z): return -mp.log(1-z)
def li2(z): return mp.polylog(2, z)
def li3(z): return mp.polylog(3, z)
def li4(z): return mp.polylog(4, z)

def G_phys(word, z=1, a=None):
    # G_phys(a_1,...,a_n; z) = int_0^z dt1/(t1-a1) int_0^{t1} dt2/(t2-a2) ...
    # recursive numerical integration; poles at 0 handled since inner G vanishes at 0
    # a: analytic continuation angles (ignored here; assume real path, letters not on (0,z))
    w = list(word)
    if len(w) == 0:
        return mp.mpf(1)
    a1 = w[0]
    f = lambda t: G_phys(w[1:], t) / (t - a1) if t != 0 else mp.mpf(0)
    # split interval to help quad
    pts = [0, z]
    return mp.quad(f, pts)

x, y = mp.mpf('0.5'), mp.mpf('0.7')
vals = {}
vals['Li22_convB_series'] = li_convB(2,2,x,y)
vals['Li31_convB_series'] = li_convB(3,1,x,y)
vals['Li21_convB_series'] = li_convB(2,1,x,y)
vals['Li12_convB_series'] = li_convB(1,2,x,y)
vals['Li11_convB_series'] = li_convB(1,1,x,y)

# candidate words (physics ordering): convB Li_{a,b}(x,y) = convA Li_{b,a}(y,x)
# convA Li_{n1,n2}(x1,x2) = (-1)^2 G_phys(0^{n2-1}, b2, 0^{n1-1}, b1; 1), b2=1/x2, b1=1/(x1*x2)
# Li22 convB: n1=n2=2, x1=y, x2=x -> word (0, 1/x, 0, 1/(x*y))
vals['G_word_Li22'] = G_phys([0, 1/x, 0, 1/(x*y)])
# Li31 convB: a=3,b=1 -> convA Li_{1,3}(y,x): n1=1,n2=3 -> word (0,0,1/x, 1/(x*y))
vals['G_word_Li31'] = G_phys([0, 0, 1/x, 1/(x*y)])
# Li21 convB: convA Li_{1,2}(y,x): n1=1,n2=2 -> word (0,1/x, 1/(x*y))
vals['G_word_Li21'] = G_phys([0, 1/x, 1/(x*y)])
# Li12 convB: convA Li_{2,1}(y,x): n1=2,n2=1 -> word (1/x, 0, 1/(x*y))
vals['G_word_Li12'] = G_phys([1/x, 0, 1/(x*y)])
# Li11 convB: convA Li_{1,1}(y,x): word (1/x, 1/(x*y))
vals['G_word_Li11'] = G_phys([1/x, 1/(x*y)])

for k, v in vals.items():
    print(k, "=", mp.nstr(v, 30))
print()
print("diff Li22:", mp.nstr(vals['Li22_convB_series'] - vals['G_word_Li22'], 8))
print("diff Li31:", mp.nstr(vals['Li31_convB_series'] - vals['G_word_Li31'], 8))
print("diff Li21:", mp.nstr(vals['Li21_convB_series'] - vals['G_word_Li21'], 8))
print("diff Li12:", mp.nstr(vals['Li12_convB_series'] - vals['G_word_Li12'], 8))
print("diff Li11:", mp.nstr(vals['Li11_convB_series'] - vals['G_word_Li11'], 8))
'''
open('steps/step_01.py','w').write(script)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','steps/step_01.py'],capture_output=True,text=True)
print(r.stdout); print(r.stderr)
open('progress.log','a').write("step_01 verify G-word vs convB series\n")
