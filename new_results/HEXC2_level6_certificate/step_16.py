for e,w,v in PSI: print('   c^%d H%s  %s'%(e,list(w),v))
print()
# Check the HEXC1 word_mpl_conv.json and w4_hpl_reductions
base = '<project>/research/08_theory'
print(open(base+'/HEXC1_hexagon_certificate/word_mpl_conv.json').read()[:1200])
