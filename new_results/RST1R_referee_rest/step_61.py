open('progress.log','a').write("step_12: independent end-to-end assembly at k=99: exact T1 (101 mpmath F-quads) = 2.18e-12 < proof bound 1.32e-9 (ratio 0.00165); every group exact<bound; T2 exact/bound 0.00039; TSV complete [0,400]. PASS\n")
import subprocess
r=subprocess.run(['grep','-n','-i','verdict\|(12)\|(13)\|(18)\|(20)\|PROVED','<project>/research/08_theory/CX6R_review/result.md'],capture_output=True,text=True)
print(r.stdout[:2500])
