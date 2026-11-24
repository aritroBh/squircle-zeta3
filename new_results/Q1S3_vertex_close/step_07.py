b3 = "<project>/research/07_push/S3_superellipsoid/claude/"
src = open(b3+"result.md").read()
# print sections 6-8
i6 = src.find("## 6")
if i6<0: i6 = src.find("# 6")
print(src[i6-200:i6+12000])
