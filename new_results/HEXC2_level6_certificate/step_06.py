base = '<project>/research/08_theory'
src = open(f'{base}/HEXC1_hexagon_certificate/cmzv6.py').read()
print(src[:2600])
print('......')
src2 = open(f'{base}/HEXC1_hexagon_certificate/hexc1_common.py').read()
print(src2[2600:4438])
