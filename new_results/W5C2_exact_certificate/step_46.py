import subprocess
print(subprocess.run(['grep','-n','def words\|def primitive\|def mul','../A3c_claude_elementary/Q/relations.py'], capture_output=True, text=True).stdout)
print(subprocess.run(['sed','-n','1,40p','../A3c_claude_elementary/Q/relations.py'], capture_output=True, text=True).stdout)