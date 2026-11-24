import sympy as sp, signal
u = sp.symbols('u', positive=True)
def handler(signum, frame): raise TimeoutError
signal.signal(signal.SIGALRM, handler)
for m in [1,2]:
    signal.alarm(90)
    try:
        r = sp.integrate(sp.sin(2*m*u)*sp.log(sp.cos(u)), u)
        print(m, "->", sp.simplify(r))
    except TimeoutError:
        print(m, "TIMEOUT")
    except Exception as e:
        print(m, "ERR", str(e)[:80])
    signal.alarm(0)
