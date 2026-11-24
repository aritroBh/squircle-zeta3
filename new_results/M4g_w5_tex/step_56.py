import subprocess, re
frag=r'''We come to $M_4$. Splitting by weight, $M_4=W_5+W_4''+W_3^{(4)}+T$, and the last three are already in closed form (proved): $W_4''$ in tetralogarithms, $W_3^{(4)}=W_3'-\frac18\zeta(3)-\frac1{24}\pi^2\lambda-\frac1{12}\lambda^3$, and $T=(\pi^2+6\lambda^2)/24$. What remains is $W_5$, the part of weight five.
$W_5$ is not classical. With $X=\lambda$ and $Y=L$, the antisymmetrised $(3,2)$ part of its symbol pairs to $84$ against $[X,[X,Y]]\otimes[X,Y]$ and to $-38\sqrt2$ against $[Y,[X,Y]]\otimes[X,Y]$. Every $\operatorname{Li}_5$ and every product of lower-weight constants pairs to zero, so no combination of $\operatorname{Li}_5$'s and products has the symbol of $W_5$.
The $(3,2)$ cobracket fixes a depth-two correction exactly. With $\operatorname{Li}_{s,t}(x,y)=\sum_{n>m>0}x^ny^m/(n^sm^t)$,
\begin{multline*}
N_5=\Bigl(779-\frac{695\sqrt2}{56}\Bigr)\operatorname{Li}_{3,2}(b,1)-672\,\operatorname{Li}_{3,2}\Bigl(\frac1{\sqrt2},1\Bigr)\\
-(374-76\sqrt2)\operatorname{Li}_{3,2}(-b,1)-\Bigl(405-\frac{4951\sqrt2}{56}\Bigr)\operatorname{Li}_{3,2}(b,-1).
\end{multline*}
The difference $W_5-N_5$ is a combination, with coefficients $p+q\sqrt2$ ($p,q$ rational), of twelve $\operatorname{Li}_5$'s, $\zeta(5)$, $L(\chi_8,5)$ and products of lower-weight constants. The twelve arguments are $b$, $-b$, $-b^2$, $1/\sqrt2$, $\frac12$, $1-1/\sqrt2$, $2-\sqrt2$, $2b$, $(2-\sqrt2)/4$, $(2+\sqrt2)/4$, $-b/2$ and $(4-3\sqrt2)/8$.
The identity holds numerically to 1596 digits, with exact agreement of the symbol and the $(3,2)$ cobracket. The coefficient of $L(\chi_8,5)$ depends on the representative, because $L(\chi_8,5)$ lies in the $\operatorname{Li}_5$ span. The full table is in the Supplement. $M_4$ is the first squircle coefficient that needs a depth-two polylogarithm.
'''
open('paper_fragment.tex','w').write(frag)
wrap='\\documentclass{article}\n\\usepackage{amsmath}\n\\begin{document}\n\\input{paper_fragment}\n\\end{document}\n'
open('fragment_test.tex','w').write(wrap)
r=subprocess.run(['nice','-n','10','/opt/homebrew/bin/tectonic','-X','compile','fragment_test.tex'],capture_output=True,text=True,timeout=600)
print('--- tectonic fragment_test ---'); print(r.stdout[-500:]); print(r.stderr[-500:])
def pages(pdf):
    data=open(pdf,'rb').read()
    m=re.findall(rb'/Count\s+(\d+)',data)
    return max(int(x) for x in m) if m else '?'
print('supplement_W5.pdf pages:', pages('supplement_W5.pdf'))
print('fragment_test.pdf pages:', pages('fragment_test.pdf'))
