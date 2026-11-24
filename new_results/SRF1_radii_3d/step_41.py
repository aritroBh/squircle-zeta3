
s3 = '<project>/research/07_push/S3_superellipsoid/claude'
r = open(s3+'/result.md').read()
print(len(r))
print(r[:7000])
