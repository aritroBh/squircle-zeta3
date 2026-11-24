import subprocess
print(subprocess.run(['bash','-c','grep -n "regular\|2m\|h = \|L = \|edge length\|circumradius" "<project>/research/08_theory/E5_edge_law/steps/step_05.py" 2>/dev/null | head -30; echo ====; sed -n "/regular/,$p" "<project>/research/08_theory/E5_edge_law/steps/step_05.py" 2>/dev/null | head -60'], capture_output=True, text=True).stdout)
