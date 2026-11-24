import os, time
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W)
with open('progress.log','a') as f:
    f.write(time.strftime('%H:%M:%S ')+"STEP3a trace1 (numpy V10 guesses, a1 in [1e-5,1e5], 121 vals): at EVERY Gamma pt QA<0, QC>0. BUT: large-a1 (a1>~1.5e4) numpy root-failure MISSED branches (mpf check at 38311,1e5 found Gamma pts with both Q11 signs + QA<0,QC>0). Re-tracing with full mpf V10 solver.\n")
    f.write(time.strftime('%H:%M:%S ')+"STEP3b launched trace2 (mpf V10 solver, a1 in [1e-7,1e6], 53 vals, dense sign-change scan + Newton). Verifies large/small-a1 ends.\n")
print("logged")
time.sleep(120)
print(open('gamma_trace2.log').read() if os.path.exists('gamma_trace2.log') else "(no log)")
