import os
print("control log:")
print(open('ctrl_z3z35.log').read() if os.path.exists('ctrl_z3z35.log') else "not yet")
print("control out tail:")
print(open('run_ctrl_z3z35.out').read() if os.path.exists('run_ctrl_z3z35.out') else "none")
