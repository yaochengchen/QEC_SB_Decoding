from pathlib import Path
import subprocess,re
import pymupdf as fitz
R=Path(__file__).resolve().parent
md=(R/'SUPPLEMENTAL_MATERIAL.md').read_text(encoding='utf-8')
md=md.replace('## S2.',r'\newpage'+'\n\n## S2.').replace(r'\(', '$').replace(r'\)', '$')
md=re.sub(r'!\[[^\]]*\]\(([^)]+)\)',lambda m:r'\begin{center}\includegraphics[width=\textwidth,height=0.55\textheight,keepaspectratio]{'+m[1]+r'}\end{center}',md)
(R/'layout.md').write_text(md,encoding='utf-8')
subprocess.run(['pandoc','layout.md','-f','markdown+raw_tex+tex_math_dollars','-s','-o','SUPPLEMENTAL_MATERIAL.tex','-V','geometry:margin=20mm','-V','fontsize=10pt','-V','mainfont=Times New Roman','--pdf-engine=xelatex'],cwd=R,check=True)
p=R/'SUPPLEMENTAL_MATERIAL.tex';t=p.read_text(encoding='utf-8').replace(r'\begin{document}',r'\usepackage{graphicx}'+'\n'+r'\begin{document}');p.write_text(t,encoding='utf-8')
for _ in range(2):
 result=subprocess.run(['xelatex','-interaction=nonstopmode','-halt-on-error','-jobname=front','SUPPLEMENTAL_MATERIAL.tex'],cwd=R,capture_output=True);assert result.returncode==0,result.stdout[-2000:]
doc=fitz.open(R/'front.pdf');assert len(doc)==2
for pv in ['0.08','0.10','0.14','0.16','0.18']:doc.insert_pdf(fitz.open(R/'figures'/f'resources_p{pv}.pdf'))
doc.save(R/'SUPPLEMENTAL_MATERIAL.pdf')
assert len(doc)==7
print('Built self-contained 7-page supplement')
