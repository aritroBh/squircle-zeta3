from towlib import *
# parent values of child (2,1): parents (3,1) [DL1], (2,2) [DL2+fold]
# level 4: both parents initial. level 5: (3,1) stable?, (2,2) stable.
print('(3,1) u(k,d):')
for dlt in range(0,4):
    print(' d=%d:'%dlt, [str(u(3,1,k,dlt)) for k in (4,5,6)])
print('(2,2) u(k,d):')
for dlt in range(0,6):
    print(' d=%d:'%dlt, [str(u(2,2,k,dlt)) for k in (4,5,6)])
# child (2,1) at level 4 vs 5
print('(2,1) u(4,d),u(5,d):', [(str(u(2,1,4,d)),str(u(2,1,5,d))) for d in range(3)])
