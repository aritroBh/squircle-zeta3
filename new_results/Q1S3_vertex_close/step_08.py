b3 = "<project>/research/07_push/S3_superellipsoid/claude/"
src = open(b3+"result.md").read()
i = src.find("## 5")
print(src[i:i+9000])
