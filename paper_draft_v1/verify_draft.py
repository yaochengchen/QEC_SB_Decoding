"""Check manuscript links, selected numerical claims, and preservation of the original paper."""
import csv
import hashlib
import json
import re
from verify_consistency import check_all
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

def rows(path):
    with (ROOT / path).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

text = (HERE / 'MANUSCRIPT.md').read_text(encoding='utf-8')
links = re.findall(r'\]\(([^)]+)\)', text)
local = [x for x in links if not x.startswith('https://')]
assert all((HERE / x).is_file() for x in local), local
original = json.loads((HERE / 'original_paper_sha256.json').read_text(encoding='utf-8-sig'))
for item in original:
    assert hashlib.sha256(Path(item['Path']).read_bytes()).hexdigest().upper() == item['Hash']

grid = rows('results/paper1_grid/grid_ler.csv')
main = [r for r in grid if r['is_main'] == 'True']
assert len(main) == 54
assert sum(int(r['shots']) for r in main) == 7_120_000
paired = rows('results/paper1_grid/grid_paired.csv')
mw = [r for r in paired if r['a'] == 'CP_local_SOE_r2' and r['b'] == 'MWPM']
assert len(mw) == 54 and all(float(r['holm_p']) < .05 for r in mw)
assert (sum(int(r['a_right_b_wrong']) for r in mw), sum(int(r['a_wrong_b_right']) for r in mw)) == (55081,17007)
bp = rows('results/paper1_extra/bposd_corr_paired.csv')
bp = [r for r in bp if r['a'] == 'CP_local_SOE_r2' and r['b'] == 'BP_OSD_corr']
assert len(bp) == 54 and all(float(r['holm_p']) < .05 for r in bp)
assert (sum(int(r['a_right_b_wrong']) for r in bp), sum(int(r['a_wrong_b_right']) for r in bp)) == (41739,21099)
budget = rows('results/paper1_extra/budget_ler.csv')
sample = [r for r in budget if int(r['d']) == 13 and float(r['p']) == .12 and r['config'] in ('sbR32S400','saC32W200')]
expected = {('sbR32S400','MAP'):1.382, ('sbR32S400','r1'):1.115, ('sbR32S400','r2'):1.096, ('saC32W200','MAP'):1.376, ('saC32W200','r1'):1.099, ('saC32W200','r2'):1.093}
assert len(sample) == 6
for r in sample:
    assert round(float(r['LER_over_REF']),3) == expected[(r['config'],r['inference'])]
matched = rows('results/paper1_extra/budget_paired_matched.csv')
comparisons = []
for a,b,n,w,l,sig in [('sbR8S400_r2','sbR32S100_r2',15,1935,1643,1),('sbR8S400_r2','saC8W100_r2',18,3169,3048,0),('sbR32S400_r2','saC32W200_r2',18,1134,1200,0)]:
    selected = [r for r in matched if r['a']==a and r['b']==b]
    assert len(selected)==n
    assert sum(int(r['a_right_b_wrong']) for r in selected)==w
    assert sum(int(r['a_wrong_b_right']) for r in selected)==l
    assert sum(float(r['holm_p'])<.05 for r in selected)==sig
    comparisons.append(dict(a=a,b=b,cells=n,wins=w,losses=l,significant=sig))
ref_numbers = [int(x) for x in re.findall(r'^(\d+)\. ',text,re.M)]
assert ref_numbers == list(range(1,33))
body = text.split('## References')[0]
citations = re.findall(r'(?<!\[)\[(\d+(?:[–,]\d+)*)\]', body)
cited = set()
for group in citations:
    for part in group.split(','):
        if '–' in part:
            lo,hi=map(int,part.split('–'))
            cited.update(range(lo,hi+1))
        else:
            cited.add(int(part))
assert cited == set(ref_numbers)
registry=(HERE/'REFERENCES_MANAGEMENT_ZH.md').read_text(encoding='utf-8')
assert [int(x) for x in re.findall(r'^### \[(\d+)\]',registry,re.M)]==ref_numbers
intro=text.split('## I. Introduction\n\n')[1].split('\n\n## II.')[0]
translation=(HERE/'INTRODUCTION_ZH.md').read_text(encoding='utf-8').split('## I. 引言\n\n')[1]
intro_blocks=len(intro.strip().split('\n\n'))
assert intro_blocks==len(translation.strip().split('\n\n'))
assert intro_blocks==json.loads((HERE/'introduction_sync.json').read_text(encoding='utf-8'))['paragraph_blocks']
assert re.findall(r'(?<!\[)\[\d+(?:[–,]\d+)*\]',intro)==re.findall(r'(?<!\[)\[\d+(?:[–,]\d+)*\]',translation)
sync=json.loads((HERE/'introduction_sync.json').read_text(encoding='utf-8'))
assert hashlib.sha256(intro.encode()).hexdigest()==sync['english_introduction_sha256']
assert text.count('\\[')==text.count('\\]')
report = dict(original_files_unchanged=len(original), local_links_valid=len(local), main_cells=len(main), main_shots=7_120_000, paired_comparisons=comparisons, references=32, introduction_translation_blocks=intro_blocks, citation_registry_consistent=True, words=len(text.split()), scope='Selected numerical and structural checks; not an independent rerun of decoding or a complete bibliography audit.')
report['extended_consistency'] = check_all(HERE)
(HERE / 'verification.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
