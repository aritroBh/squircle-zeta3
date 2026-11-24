import subprocess
print(subprocess.run(['bash','-c','grep -rl "656.4153619676247544766476" "<project>/research/" 2>/dev/null | head; grep -rn "656.4153619676" "<project>/research/08_theory/D17"* 2>/dev/null | head -5'], capture_output=True, text=True).stdout)
