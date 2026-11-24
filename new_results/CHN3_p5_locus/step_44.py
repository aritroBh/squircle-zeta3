import subprocess, time, os
# 1. Log progress
with open('progress.log','a') as f:
    f.write("CHN3 step1: V10={C4=Q7'=0}. Res_a5 (a3=1; lcs 4,64 const) = 16*a1*(a1+a3)*G, G irred deg22 (2138 terms). "
            "Positive part: only via {G=0} (a1, a1+a3 factors positive-infeasible). V10 birational to {G=0} via a5=-B/A "
            "(deg-1 subresultant). EXACT certified positive point (1,1,1,a4*,a5*=-B/A), a4*=root G(1,1,.) in (1,2), "
            "a5* in (0.07051999939180706,...) >0 via exact interval arithmetic; second point a4 in (5,6), a5~6.5304. "
            "Files: R_a1a2a4.txt, R_factors.txt, subres1_A/B.txt, subres0.txt\n")
print("logged")

# 2. Launch engine W<=12 build (heavy, nice)
PY = '<project>/.venv/bin/python'
logf = open('engine12.log','w')
p = subprocess.Popen(['nice','-n','10',PY,'run_engine12.py'], stdout=logf, stderr=subprocess.STDOUT,
                     cwd='<project>/research/08_theory/CHN3_p5_locus')
print("engine12 PID:", p.pid)
time.sleep(30)
print(open('engine12.log').read())
