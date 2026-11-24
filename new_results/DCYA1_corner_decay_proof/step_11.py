import subprocess
print("=== corner_series.py ===")
print(subprocess.run(['cat', '<project>/research/10_frontier/G1_zeta3_any_angle/corner_series.py'], capture_output=True, text=True).stdout)
