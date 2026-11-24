paper = open('../../../paper/ieee/main.tex').read()
print(len(paper))
# Find the pure weight subsection
i = paper.find('Pure Weight')
print(paper[i-2000:i+8000])
