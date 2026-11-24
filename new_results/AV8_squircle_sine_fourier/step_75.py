import sympy as sp, signal, time
y = sp.symbols('y', positive=True)
s2=y**2; c2=1-y**2; ls=sp.log(y); lc2=sp.log(1-y**2)
A = 2*s2*ls - (s2-1)*lc2
H2 = (8*s2**2*ls**2 - 8*s2*(s2-1)*(ls-1)*ls - 4*c2**2*lc2 + A**2)/(16*c2**sp.Rational(5,2))
# split into pieces and integrate each with time limit
pieces = sp.expand(H2).args
print("num pieces:", len(pieces))
def handler(signum, frame): raise TimeoutError
signal.signal(signal.SIGALRM, handler)
results = []
t0=time.time()
for i,pc in enumerate(pieces):
    signal.alarm(60)
    try:
        r = sp.integrate(pc, y)
        results.append((pc, r, "OK"))
    except TimeoutError:
        results.append((pc, None, "TIMEOUT"))
    except Exception as e:
        results.append((pc, None, str(e)[:50]))
    signal.alarm(0)
ok = sum(1 for r in results if r[2]=="OK")
print(f"integrated {ok}/{len(pieces)} in {time.time()-t0:.0f}s")
for pc,r,st in results:
    print(st, "|", pc, "->", str(r)[:120])
