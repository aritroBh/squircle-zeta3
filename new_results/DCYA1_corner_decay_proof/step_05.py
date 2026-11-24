import subprocess
# Read the corner decay law section (1455-1485)
print("=== Corner decay law section (1455-1485) ===")
print(subprocess.run(['sed', '-n', '1455,1485p', '<project>/paper/ieee/main.tex'], capture_output=True, text=True).stdout)
print("\n=== M_k(alpha) definition section (1377-1410) ===")
print(subprocess.run(['sed', '-n', '1377,1410p', '<project>/paper/ieee/main.tex'], capture_output=True, text=True).stdout)
