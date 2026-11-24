import subprocess
print(subprocess.run(['bash','-c','ls "<project>/paper/ieee/" 2>/dev/null | head; grep -n "thm:pip" "<project>/paper/ieee/main.tex" | head'], capture_output=True, text=True).stdout)
