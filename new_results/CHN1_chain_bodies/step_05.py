import os
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
md='<project>/research/02_ideas/lenses/experimental.md'
lines=open(md).read().splitlines()
print('\n'.join(lines[124:215]))
