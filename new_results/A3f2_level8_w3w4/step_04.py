import os
Q='<project>/research/08_theory/A3c_claude_elementary/Q'
print(os.listdir(Q))
print(open(Q+'/progress.log').read()[:4000])
