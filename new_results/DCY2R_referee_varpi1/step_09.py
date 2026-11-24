import subprocess, os
# progress.log and result.md scaffolding in my cwd
open('progress.log','w').write("DCY2R referee log\n")
open('result.md','w').write("# DCY2R referee report (in progress)\n")
# find the paper
print(subprocess.run(['bash','-c','ls ../../../paper/ieee/ 2>/dev/null; ls ../../../paper 2>/dev/null'],capture_output=True,text=True).stdout)
print(subprocess.run(['bash','-c','grep -n "varpi" ../../../paper/ieee/main.tex | head -50'],capture_output=True,text=True).stdout)
