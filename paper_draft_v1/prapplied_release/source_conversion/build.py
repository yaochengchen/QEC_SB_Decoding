"""Faithful Markdown to LaTeX release; run with Python, Pandoc and XeLaTeX."""
from pathlib import Path
import re, shutil, subprocess, hashlib, json
ROOT = Path(__file__).resolve().parent
DRAFT = ROOT.parent.parent
md = (DRAFT/'MANUSCRIPT.md').read_text(encoding='utf-8')
source = md.replace(r'\(', '$').replace(r'\)', '$').replace(r'\[', '$$').replace(r'\]', '$$')
(ROOT/'figures').mkdir(exist_ok=True)
images=[]
def image(m):
    origin=(DRAFT/m[2]).resolve(); target=ROOT/'figures'/origin.name
    shutil.copy2(origin,target)
    images.append({'source':str(origin),'copy':str(target),'sha256':hashlib.sha256(origin.read_bytes()).hexdigest()})
    return f'![{m[1]}](figures/{origin.name})'
source=re.sub(r'!\[([^\]]*)\]\(([^)]+)\)', image, source)
# Local evidence links are one directory deeper in this release. No prose changes.
source=re.sub(r'(?<!!)\[([^\]]+)\]\(([^)]+)\)',lambda m: m[0] if re.match(r'(?:https?://|mailto:)',m[2]) else f'[{m[1]}](../{m[2]})',source)
(ROOT/'source.md').write_text(source,encoding='utf-8')
(ROOT/'image_manifest.json').write_text(json.dumps(images,indent=2),encoding='utf-8')
subprocess.run(['pandoc','source.md','-f','markdown-implicit_figures+tex_math_dollars','-s','-t','latex','-o','manuscript.tex','--pdf-engine=xelatex','-V','geometry:margin=22mm','-V','fontsize=10pt','-V','mainfont=Times New Roman','-H','header.tex'],cwd=ROOT,check=True)
tex=(ROOT/'manuscript.tex').read_text(encoding='utf-8')
# Keep each original figure and its complete manually numbered caption together.
pattern=r'(\\pandocbounded\{\\includegraphics[^\n]+\})\s*\n\s*(\\textbf\{Figure (S?\d+)\.[\s\S]*?)(?=\n\n)'
def fig(m):
    body=m[1].replace('keepaspectratio,','width=\\linewidth,height=0.62\\textheight,keepaspectratio,')
    return '\\begin{figure}[H]\n\\centering\n'+body+'\n\\par\\smallskip\n\\begin{minipage}{\\linewidth}\n'+m[2]+'\n\\end{minipage}\n\\label{fig:'+m[3]+'}\n\\end{figure}'
tex,n=re.subn(pattern,fig,tex); assert n==5,n
# Title and heading levels only; original wording and explicit section numbers retained.
tex=re.sub(r'\\section\{([\s\S]*?)\}\\label\{[^}]+\}',lambda m:'\\begin{center}\n{\\LARGE\\bfseries '+m[1]+'\\par}\n\\end{center}',tex,count=1)
tex=tex.replace(r'\subsection{',r'\section{').replace(r'\subsubsection{',r'\subsection{')
# Center the supplied author and affiliation block beneath the title.
tex, author_blocks = re.subn(r'(\\textbf\{Jiapei Zhuang[\s\S]*?)(?=\\section\{Abstract\})', lambda m: '\\begin{center}\n\\small\n'+m[1]+'\\end{center}\n\n', tex)
assert author_blocks == 1, author_blocks
# Break long math lists in the evaluation table and tune the structural table.
def table(m):
    t=m[0]
    if 'Main evaluation grid' in t:
        t=t.replace(',',r',\allowbreak ')
    if 'Maximum mechanism incidence' in t:
        t=t.replace('0.1111','0.1600').replace('0.1481','0.1400')
        t='{\\small\\setlength{\\tabcolsep}{4pt}\n'+t+'\n}'
    return t
tex=re.sub(r'\\begin\{longtable\}[\s\S]*?\\end\{longtable\}',table,tex)
# Keep each numbered table and its full caption on one page.
tex=re.sub(r'(\\textbf\{Table (?:1|2|C1)\.[\s\S]*?\\end\{longtable\}(?:\n\})?)', lambda m: '\\par\n\\noindent\\begin{minipage}{\\linewidth}\n'+m[1]+'\n\\end{minipage}\n\\par', tex)
# Explicit existing numbers are used as cross-reference labels, without renumbering.
tex=re.sub(r'\\label\{fig:(S?\d+)\}',lambda m:r'\fixedlabel{fig:'+m[1]+'}{'+m[1]+'}',tex)
tex=re.sub(r'\\tag\{(A?\d+)\}',lambda m:m[0]+r'\label{eq:'+m[1]+'}',tex)
tex=re.sub(r'(?=\\textbf\{Table (1|2|C1)\.)',lambda m:r'\fixedlabel{tab:'+m[1]+'}{'+m[1]+'}',tex)
def crossref(m):
    kind,num=m[1],m[2]
    prefix='eq' if kind=='Eq.' else ('tab' if kind=='Table' else 'fig')
    number=num.strip('()')
    if prefix=='fig' and number not in ['1','2','3','4','5']: return m[0]
    return r'\hyperref['+prefix+':'+number+']{'+m[0]+'}'
tex=re.sub(r'(?<!\{)\b(Fig\.|Figure|Table|Eq\.)\s+(\(A?\d+\)|S?\d+|C1)\b',crossref,tex)
# Give references PDF destinations and link existing numeric citations, preserving numbering.
start=tex.index(r'\section{References}')
body,bib=tex[:start],tex[start:]
count=[0]
def refitem(m):
    count[0]+=1
    return m[0]+'\\hypertarget{ref:'+str(count[0])+'}{} '
bib=re.sub(r'\\item\s*',refitem,bib)
assert count[0]==32,count
# Pandoc escapes square brackets with {[} and {]}.
def cite(m):
    nums=m[1].split(',')
    return '{[}'+','.join('\\hyperlink{ref:'+v.strip()+'}{'+v+'}' for v in nums)+'{]}'
body=re.sub(r'\{\[\}([0-9]+(?:,[0-9]+)*)\{\]\}',cite,body)
bib=bib.replace(r'\def\labelenumi{\arabic{enumi}.}', r'\def\labelenumi{[\arabic{enumi}]}')
tex=body+bib
(ROOT/'manuscript.tex').write_text(tex,encoding='utf-8')
