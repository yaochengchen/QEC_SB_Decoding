from pathlib import Path
import re,shutil,subprocess,json
R=Path(__file__).resolve().parent
subprocess.run(['python',str(R/'source_conversion/build.py')],check=True)
source=(R/'source_conversion/manuscript.tex').read_text(encoding='utf-8')
shutil.copytree(R/'source_conversion/figures',R/'figures',dirs_exist_ok=True)
body=source.split(r'\section{I. Introduction}',1)[1].split(r'\section{References}',1)[0]
body=r'\section{Introduction}'+body
abstract=source.split(r'\section{Abstract}\label{abstract}',1)[1].split(r'\section{I. Introduction}',1)[0].strip()
# Native section numbering, preserving the visible manuscript numbering.
body=re.sub(r'\\section\{(?:II|III|IV|V|VI)\. ',r'\\section{',body)
body=re.sub(r'\\subsection\{[A-F]\. ',r'\\subsection{',body)
body=body.replace(r'\section{Appendix A.',r'\appendix'+'\n'+r'\renewcommand{\thesubsection}{\thesection.\arabic{subsection}}'+'\n'+r'\section{')
body=re.sub(r'\\section\{Appendix [BC]\. ',r'\\section{',body)
body=re.sub(r'\\subsection\{[AC]\.\d+\. ',r'\\subsection{',body)
# Float bodies and caption wording are preserved; native counters retain their IDs.
def figure(m):
 s=m[0]; num=re.search(r'\\fixedlabel\{fig:([^}]+)',s)[1]
 image=re.search(r'\\includegraphics[^\n]*\{(figures/[^}]+)\}',s)[1]
 cap=re.search(r'\\textbf\{Figure '+num+r'\. ([\s\S]*?)\}([\s\S]*?)\\end\{minipage\}',s)
 assert cap,num
 content=cap[1]+cap[2]
 return '\n'+(r'\renewcommand{\thefigure}{S\arabic{figure}}\setcounter{figure}{0}'+'\n' if num=='S1' else '')+r'\begin{figure*}[t]'+'\n'+r'\centering\includegraphics[width=\textwidth,height=0.48\textheight,keepaspectratio]{'+image+'}\n'+r'\caption{'+content.strip()+'}'+r'\label{fig:'+num+'}\n'+r'\end{figure*}'+'\n'
body,n=re.subn(r'\\begin\{figure\}\[H\][\s\S]*?\\end\{figure\}',figure,body);assert n==5
# Convert Pandoc longtables to standard tabulars, retaining every cell.
def tabular(s):
 s=s.replace(r'\raggedleft',r'\raggedright')
 s=s.replace(r'\begin{longtable}[]',r'\begin{tabular}').replace(r'\end{longtable}',r'\bottomrule\end{tabular}')
 s=s.replace('\\endhead\n\\bottomrule\\noalign{}\n\\endlastfoot','')
 s=re.sub(r'>\{\\ragged(?:right|left)\\arraybackslash\}', '', s)
 s=re.sub(r'p\{\(\\linewidth - (\d+)\\tabcolsep\) \* \\real\{([0-9.]+)\}\}',lambda m:'p{'+str(float(m[2])*0.94)+r'\linewidth}',s)
 return s.replace(r'\noalign{}','')
def table(m):
 num,caption,tablebody=m[1],m[2],m[3]
 env='table' if num=='2' else 'table*'
 # Native labels use manuscript Arabic table numbering.
 return ('\n'+(r'\renewcommand{\thetable}{C\arabic{table}}\setcounter{table}{0}'+'\n' if num=='C1' else '')+r'\begin{'+env+'}[t]\n'+r'\caption{'+caption.strip()+'}'+r'\label{tab:'+num+'}\n'+r'\centering\small'+'\n'+tabular(tablebody)+'\n'+r'\end{'+env+'}\n')
pattern=r'\\par\s*\\noindent\\begin\{minipage\}\{\\linewidth\}\s*\\fixedlabel\{tab:([^}]+)\}\{[^}]+\}\\textbf\{Table [^ ]+\. ([\s\S]*?)\}\s*([\s\S]*?\\end\{longtable\}(?:\s*\})?)\s*\\end\{minipage\}\s*\\par'
# Capture prose between bold heading and table as part of native caption.
def tablewrap(m):
 rest=m[3]; pos=rest.index(r'\begin{longtable}')
 # Table C1 has a small-font grouping just before the tabular.
 prefix=rest[:pos]; prefix=re.sub(r'\{\\small\\setlength\{\\tabcolsep\}\{4pt\}\s*$','',prefix)
 tbl=rest[pos:];tbl=re.sub(r'\s*\}\s*$','',tbl) if tbl.rstrip().endswith('}\n}') else tbl
 class Match:
  def __getitem__(self,i):return [None,m[1],m[2]+(' '+prefix.strip() if prefix.strip() else ''),tbl][i]
 return table(Match())
body,n=re.subn(pattern,tablewrap,body);assert n==2,n
# The unnumbered reproducibility index remains unnumbered and spans two columns.
body=re.sub(r'\\begin\{longtable\}[\s\S]*?\\end\{longtable\}',lambda m:r'\begin{table*}[t]\centering\small'+'\n'+tabular(m[0])+'\n'+r'\end{table*}',body)
# Reflow long archive filenames without changing displayed characters.
body=re.sub(r'\\href\{([^}]+)\}\{((?:[^{}]|\\[{}_])*)\}',lambda m:r'\href{'+m[1]+'}{'+m[2].replace('/',r'/\allowbreak ').replace(r'\_',r'\_\allowbreak ')+'}',body)
# Restore normal native citation handling.
body=re.sub(r'\{\[\}((?:\\hyperlink\{ref:\d+\}\{\d+\},?)+)\{\]\}',lambda m:r'\cite{'+','.join('ref'+n for n in re.findall(r'\\hyperlink\{ref:(\d+)\}',m[1]))+'}',body)
bib=source.split(r'\begin{enumerate}',1)[1].split(r'\end{enumerate}',1)[0]
items=re.findall(r'\\item\s*\\hypertarget\{ref:(\d+)\}\{\}\s*([\s\S]*?)(?=\\item|\Z)',bib);assert len(items)==32
bib='\n'.join(r'\bibitem{ref'+n+'} '+v.strip() for n,v in items)
preamble=r"""\documentclass[aps,prapplied,reprint,superscriptaddress,amsmath,amssymb,floatfix,nofootinbib]{revtex4-2}
\usepackage{fontspec}
\usepackage{unicode-math}
\setmainfont{Times New Roman}
\setmathfont{Latin Modern Math}
\usepackage{graphicx,booktabs,calc}
\usepackage{xurl}
\usepackage[hidelinks,unicode]{hyperref}
\urlstyle{same}
\newcommand{\tightlist}{\raggedright\setlength{\itemsep}{0pt}\setlength{\parskip}{0pt}}
\renewcommand{\thetable}{\arabic{table}}
\setlength{\emergencystretch}{2em}
\begin{document}
\title{Quantum-Inspired Search Meets Quantum Error Decoding: Constraint-Preserving Simulated Bifurcation with Local Degeneracy Aggregation}
\author{Jiapei Zhuang}
\affiliation{Quantum Science Center of Guangdong-Hong Kong-Macao Greater Bay Area, Shenzhen 518045, China}
\author{Yaocheng Chen}
\email{yaochengchen@usp.br}
\affiliation{Instituto de Astronomia, Geofísica e Ciências Atmosféricas (IAG-USP), Universidade de São Paulo, Rua do Matão 1226, São Paulo, SP 05508-090, Brazil}
"""
# Line breaks only: no symbol or mathematical term is changed.
eq_index=[0]
def reflow(m):
 i=eq_index[0];eq_index[0]+=1;x=m[1]
 if i not in [2,3,5,7,11,14,17,19]:return m[0]
 tag=re.search(r'\\tag\{[^}]+\}\\label\{[^}]+\}',x)
 tail=tag[0] if tag else ''
 if tag:x=x.replace(tail,'')
 if i==2:x=x.replace(r',\quad',',\\\\',1)
 elif i in [3,5,7,11,14,17]:x=x.replace(r'\qquad',r'\\')
 elif i==19:x=x.replace('+\\sum',r'\\'+'+\\sum',1)
 return '\\[\n\\begin{gathered}\n'+x.strip()+'\n\\end{gathered}'+('\n'+tail if tail else '')+'\n\\]'
body=re.sub(r'\\\[([\s\S]*?)\\\]',reflow,body)

tex=preamble+r'\begin{abstract}'+'\n'+abstract+'\n'+r'\end{abstract}\maketitle'+'\n'+body+'\n'+r'\begin{thebibliography}{32}'+'\n'+bib+'\n'+r'\end{thebibliography}\end{document}'
(R/'manuscript.tex').write_text(tex,encoding='utf-8')
for i in range(3):
 run=subprocess.run(['xelatex','-interaction=nonstopmode','-halt-on-error','manuscript.tex'],cwd=R,capture_output=True,text=True,encoding='utf-8',errors='replace')
 (R/f'build_pass_{i+1}.txt').write_text(run.stdout,encoding='utf-8')
 if run.returncode:raise RuntimeError(run.stdout[-3500:])
print('REVTeX PDF built.')
