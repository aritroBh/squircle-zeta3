import subprocess
print(subprocess.run(['bash','-c','grep -n "thm:hexagon5\|20187\|30720\|2835" "<project>/paper/ieee/main.tex" | head -20'], capture_output=True, text=True).stdout)
