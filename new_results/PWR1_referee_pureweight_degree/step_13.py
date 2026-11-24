paper = open('../../../paper/ieee/main.tex').read()
i = paper.find('Fourier Series Meet the Squircle')
print(paper[i:i+6500])
