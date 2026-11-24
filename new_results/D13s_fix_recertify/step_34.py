
import os, subprocess, time
WORK = 'research/08_theory/D13s_fix_recertify/D13s_work'
FIX = os.path.join(WORK, 'fixed')
PY='<home>)/.venv/bin/python'
with open(os.path.join(WORK,'progress.log'),'a') as f:
    f.write("D13s turn2-6: B2-fix math PROVED (f(x)=1-atan(x)/x nondecr: f'=g/x^2, g(0)=0, g'=2x^2/(1+x^2)^2>=0; sympy printed).\n"
            " H2 PROVED: x^2|g(x)|<=1 via sHc(s)=1-sum c_k s/(s+2k) + 2Re w_k>=|w_k| (3x^2+6rt2 k x+4k^2>0); AM-GM 1+rt2x+x^2>=(2+rt2)x;\n"
            "   => I0 <= (2+rt2)^{-N/2} 2/(N-2); truth ratios I0/proved ~0.27 (N=100..1e4 printed).\n"
            " H3: identity RePhi0+th = -(N/2)S(t) - (pi rt2/(4t))(rho-1)^2 verified (1.8e-31); upper h3 box valid ALL t>0 (margin 3.9e-21);\n"
            "   lower box: H(t)>0 for t<=4.54 (zero at 4.5447, N=0.288); H'(t)=t^3 Hp3(t), Hp3 manifestly regular, min on [0,0.127]=0.4027.\n"
            "   NOTE: review's margin formula is NOT a valid lower bound for small t (its bracket exceeds true margin by ~(1/6)t^4(N/2));\n"
            "   review's 'breg_h-1/12 exact' claim FALSE (their step_32 prints deviation up to 0.044; step_31.py IndentationError, never ran).\n"
            " H1: derived EXACT regular form Btrue = rho*B1/(2+rt2 t), B1 = 1/2+rt2 s1+s1^2-rho^2 A1n (s1=-t/((1+rt2 t)(1+rho)),\n"
            "   A1n in [zz2/3-zz2^2/5, zz2/3]/t^2): worst rel resid 8.5e-37 over t in [1e-12,1]; sup|Btrue| final cell = 1/6.\n"
            " fixed/step_38.py written (B2+H2 fixes). A1 box containment: 14112 pts, worst viol -9.9e-7. Bh[A1'/y] exact (1.9e-53).\n"
            " z-box containment with fixed code: 0 misses / 369 sampled N over 9 cells [100,1.5e9].\n")
t0 = time.time()
r = subprocess.run([PY, 'step_18.py'], cwd=FIX, capture_output=True, text=True, timeout=3600)
print(f"({time.time()-t0:.1f}s)")
print(r.stdout); print(r.stderr[-1200:] if r.stderr else "")
