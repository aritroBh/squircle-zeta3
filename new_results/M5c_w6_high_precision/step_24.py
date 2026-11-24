import subprocess
work='<project>/research/08_theory/M5c_w6_high_precision'
r=subprocess.run(['bash','-lc',f'cat "{work}/w6_hp.out"; echo ---; ps -p $(cat "{work}/w6_hp.pid") -o pid,etime,%cpu,comm; echo ---; tail -3 "{work}/progress.log"'],capture_output=True,text=True)
print(r.stdout); print(r.stderr)
