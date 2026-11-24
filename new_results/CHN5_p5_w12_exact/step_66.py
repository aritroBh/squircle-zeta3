import os, subprocess, time
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W)
with open('progress.log','a') as f:
    f.write(time.strftime('%H:%M:%S ')+"CAMPAIGN DONE: all 924 pts evaluated (g=vw+vwr, m0), 5 chunks x ~1200s. [eval12_*.pkl]\n")
    f.write(time.strftime('%H:%M:%S ')+"STEP1 PROVED: m0 contributes EXACTLY 0 to A,B,CLASS,Z at all 924 pts. Interpolated quotients QA(186 orbit coeffs),QB(182),QC(230),QZ(115). D = f - a2a3a4(Sum a)*Q8 == 0 at ALL 924 pts for ALL FOUR coords. det12 rank 924 mod two primes => A=L*QA, B=L*QB, CLASS=L*QC, Z=L*QZ PROVED (L=a2a3a4(Sum a)), reversal-inv deg-8 quotients. [interp12.log, Q8_quotients.pkl]\n")
print("logged")
proc=subprocess.Popen(['nice','-n','10','<project>/.venv/bin/python','run_ideal.py'],
                      cwd=W, stdout=open(W+'/run_ideal.out','w'), stderr=subprocess.STDOUT)
proc.wait()
print(open(W+'/ideal.log').read())
print(open(W+'/run_ideal.out').read())
