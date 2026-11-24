
src = open("reflib.py").read()
i0 = src.index("def hpl_value(w, x):")
clean = '''def hpl_value(w, x, _memo={}):
    """Regularized HPL value H_w(x), letters 0,±1 (0 = dt/t, ±1 = dt/(1∓t)), first letter
    outermost. Trailing 0's (log x divergence) are removed by MY shuffle identity
    H_u * H_0 = sum_v mult_v H_v  (H_0(x) = log x, T = 0 convention); convergent words are
    evaluated as truncated multiple sums (geometric for |x| <= 1/2)."""
    w = tuple(w)
    key = (w, mp.mp.dps)
    if key in _memo: return _memo[key]
    if not w: return mp.mpc(1)
    L = mp.log(x)
    if w == (0,):
        return L
    if w[-1] != 0:
        letters = tuple(('T',) if j == 0 else ('A', mp.mpc(j)) for j in w)
        val = msum(_blocks(letters), mp.mpf(x))
    else:
        u = w[:-1]
        shuf = dict(sh(u, (0,)))        # H_u ш H_(0) = sum_v mult_v H_v
        m = shuf.pop(w)                 # multiplicity of w = (trailing zero run)+1
        tot = L * hpl_value(u, x)
        for v, mult in shuf.items():
            tot -= mult * hpl_value(v, x)
        val = tot / m
    _memo[key] = val
    return val
'''
src = src[:i0] + clean
open("reflib.py", "w").write(src)
import subprocess
r = subprocess.run(["nice", "-n", "10", "<project>/.venv/bin/python", "steps/step_03.py"], capture_output=True, text=True, timeout=3000)
print(r.stdout[-1800:])
print(r.stderr[-1500:] if r.returncode else "OK")
