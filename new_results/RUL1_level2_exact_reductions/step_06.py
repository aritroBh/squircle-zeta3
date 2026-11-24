import subprocess
print(subprocess.run(['bash','-c','for d in K1f_hexagon_order5 K1g_hexagon_order5_assembly E5_edge_law E5c_nonlocal_z5 E5d_nonlocal_pslq E5e_pslq_340; do echo "##### $d"; ls "<project>/research/08_theory/$d" | head -20; done'], capture_output=True, text=True).stdout)
