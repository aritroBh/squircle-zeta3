
import subprocess
r = subprocess.run(['<project>/.venv/bin/python','pslq2.py','8','400'],capture_output=True,text=True,timeout=1600)
print(r.stdout, r.stderr[-200:])
with open('progress.log','a') as f:
    f.write("step_12c/13/pslq2: PSLQ record c5(100,150 match), c6(150,200 match maxc3600), c7(300 match maxc30240; 100,200 overfit caution); incomplete-basis controls behave correctly\n")
    f.write("step_13b/fix3: top-weight identity P_n(1,1)=2P_n(1,0) verified n=1..8; deg-n part of e_n = (-1)^n(lv^n+lw^n)/(2n!)\n")
    f.write("fix5: duality V(eta)=V(1/eta) verified 50-69 digits (well-conditioned theta-form)\n")
