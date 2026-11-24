import subprocess
print(subprocess.run(['bash','-c','grep -n "J_" "<project>/research/08_theory/D12_pi_p/result.md" | head -40; echo ====; sed -n "/J-evaluations/,+6p" "<project>/research/08_theory/D12_pi_p/result.md"'], capture_output=True, text=True).stdout)
