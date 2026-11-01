from pathlib import Path
import re,json,hashlib,subprocess
import pymupdf as fitz
R=Path(__file__).resolve().parent
md=(R.parent/'MANUSCRIPT.md').read_text(encoding='utf-8')
tex=(R/'manuscript.tex').read_text(encoding='utf-8')
log=(R/'manuscript.log').read_text(encoding='utf-8')
def maths(s):
    return [re.sub(r'\s+','',re.sub(r'\\label\{[^}]+\}','',a or b).replace(r'\allowbreak','')) for a,b in re.findall(r'\\\((.*?)\\\)|\\\[(.*?)\\\]',s,re.S)]
assert maths(md)==maths(tex)
tags=re.findall(r'\\tag\{([^}]+)\}',tex)
assert tags==[str(i) for i in range(1,11)]+['A1','A2']
refs=re.findall(r'\\hypertarget\{ref:(\d+)\}',tex)
assert refs==[str(i) for i in range(1,33)]
raw=subprocess.run(['pandoc','source.md','-f','markdown-implicit_figures+tex_math_dollars','-s','-t','latex'],cwd=R,capture_output=True,text=True,encoding='utf-8',check=True).stdout
rawbib=raw[raw.index(r'\begin{enumerate}',raw.index('{References}')):]
finalbib=tex[tex.index(r'\begin{enumerate}',tex.index('{References}')):]
finalbib=re.sub(r'\\hypertarget\{ref:\d+\}\{\} ','',finalbib)
rawbib=rawbib.replace(r'\def\labelenumi{\arabic{enumi}.}', r'\def\labelenumi{[\arabic{enumi}]}')
assert rawbib==finalbib
assert not re.search(r'[\u4e00-\u9fff]', md.split('## References',1)[1])
manifest=json.loads((R/'image_manifest.json').read_text())
assert len(manifest)==5
for item in manifest:
    assert hashlib.sha256(Path(item['source']).read_bytes()).hexdigest()==item['sha256']
    assert hashlib.sha256(Path(item['copy']).read_bytes()).hexdigest()==item['sha256']
labels=set(re.findall(r'\\(?:label|fixedlabel)\{([^}]+)\}',tex))
links=re.findall(r'\\hyperref\[([^\]]+)\]',tex)
assert set(links)<=labels
citations=re.findall(r'\\hyperlink\{ref:(\d+)\}',tex)
assert set(citations)==set(refs)
assert not re.search(r'Overfull|undefined|Missing character|^!',log,re.M)
original=json.loads((R.parent/'original_paper_sha256.json').read_text(encoding='utf-8-sig'))
for item in original:
    assert hashlib.sha256(Path(item['Path']).read_bytes()).hexdigest().upper()==item['Hash'].upper()
d=fitz.open(R/'manuscript.pdf')
assert sum(len(p.get_images()) for p in d)==5
for p in d:
    for b in p.get_text('dict')['blocks']:
        x0,y0,x1,y1=b['bbox'];assert x0>=0 and y0>=0 and x1<=p.rect.width+1 and y1<=p.rect.height+1
    for link in p.get_links():
        if link['kind']==fitz.LINK_GOTO:assert 0<=link['page']<len(d)
figure_pages={}
for num in ['1','2','3','4','5']:
    found=[i+1 for i,p in enumerate(d) if f'Figure {num}.' in p.get_text()]
    assert len(found)==1,(num,found)
    assert d[found[0]-1].get_images()
    figure_pages[num]=found[0]
result={'pages':len(d),'math_entries_identical':len(maths(md)),'equation_numbers':tags,'references_unchanged':len(refs),'citation_links':len(citations),'figure_table_crossref_links':len(links),'figures_unchanged':len(manifest),'figure_pages':figure_pages,'original_paper_files_unchanged':len(original),'overfull_or_unresolved_reference_warnings':0,'underfull_line_warnings':len(re.findall(r'Underfull',log)),'manuscript_sha256':hashlib.sha256((R.parent/'MANUSCRIPT.md').read_bytes()).hexdigest(),'pdf_sha256':hashlib.sha256((R/'manuscript.pdf').read_bytes()).hexdigest(),'scope':'Conversion/numbering/layout checks; no new experiments or external bibliography audit.'}
(R/'release_verification.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
