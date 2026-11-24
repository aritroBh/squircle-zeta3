import subprocess
print(subprocess.run(['bash','-c','ls "<project>/research/08_theory/A3c_claude_elementary/Q/" 2>/dev/null; echo ====; sed -n 1,120p "<project>/research/08_theory/A3c_claude_elementary/Q/cmzv.py"'], capture_output=True, text=True).stdout)
