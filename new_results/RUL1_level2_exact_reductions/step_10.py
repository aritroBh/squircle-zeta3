import subprocess
print(subprocess.run(['bash','-c','cat "<project>/research/08_theory/E5c_nonlocal_z5/result.md"; echo "########## EDGE LAW"; cat "<project>/research/08_theory/E5_edge_law/result.md"'], capture_output=True, text=True).stdout)
