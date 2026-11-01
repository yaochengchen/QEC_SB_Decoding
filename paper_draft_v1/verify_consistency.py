"""Structural consistency checks; not a semantic translation or external-source audit."""
from pathlib import Path
import csv
import hashlib
import json
import re

PAIRS = [
    ('Abstract', 'ABSTRACT_ZH.md'), ('I', 'INTRODUCTION_ZH.md'),
    ('II', 'SECTION_II_ZH.md'), ('III', 'SECTION_III_ZH.md'),
    ('IV', 'SECTION_IV_ZH.md'), ('V', 'SECTION_V_ZH.md'),
    ('VI', 'SECTION_VI_ZH.md'), ('A', 'APPENDIX_A_ZH.md'),
]

def citations(text):
    # Exclude [[n,k,d]] code parameters, Markdown links, and inline code.
    text = re.sub(r'`[^`]*`', '', text)
    return re.findall(r'(?<!\[)\[(\d+(?:[,–]\d+)*)\](?![\](])', text)

def citation_locations(text):
    result = {}
    section = ''
    location = ''
    for line in text.splitlines():
        if line == '## References':
            break
        if line.startswith('## '):
            label = line[3:]
            section = label.split('.')[0].replace('Appendix ', '')
            location = section
        elif line.startswith('### '):
            sub = line[4:].split('.')[0]
            location = line[4:].split()[0].rstrip('.') if section in ['A','B','C'] else section+' '+sub
        for group in citations(line):
            for part in group.split(','):
                ids = range(int(part.split('–')[0]), int(part.split('–')[1])+1) if '–' in part else [int(part)]
                for n in ids:
                    if location not in result.setdefault(n, []):
                        result[n].append(location)
    return result

def sections(text):
    chunks = re.split(r'^## ', text, flags=re.M)[1:]
    return {chunk.split('\n',1)[0].split('.')[0].replace('Appendix ',''): chunk for chunk in chunks}

def normalized_math(text):
    text = re.sub(r'\\(?:qquad|quad)', '', text)
    text = re.sub(r'\s+', '', text)
    # Punctuation is cosmetic, decimal points between digits are not.
    text = text.replace(',', '')
    text = re.sub(r'(?<!\d)\.|\.(?!\d)', '', text)
    return text

def assert_text_pair(english, chinese):
    get_math = lambda t: [normalized_math(v) for v in re.findall(r'\\\[(.*?)\\\]', t, re.S)]
    assert get_math(english) == get_math(chinese), 'Display mathematics differ'
    images = lambda t: re.findall(r'!\[[^\]]*\]\(([^)]+)\)', t)
    assert images(english) == images(chinese), 'Image sources differ'
    assert citations(english) == citations(chinese), 'Citation order differs'
    assert re.findall(r'\\tag\{([^}]+)\}',english) == re.findall(r'\\tag\{([^}]+)\}',chinese), 'Equation tags differ'

def source_refs(text):
    return {int(n):v.strip() for n,v in re.findall(r'^(\d+)\. (.+)$',text.split('## References')[1],re.M)}

def reference_index(text):
    refs = source_refs(text)
    lines=['## 当前引用位置索引','',
        '以下位置根据当前英文正文生成；I–VI表示正文节，C.1等表示附录小节。同一位置内重复引用合并列出。历史核验记录不作为当前位置索引。此表仅核对内部引用对应关系，不代表重新核验外部出版状态。','',
        '| 编号 | 当前引用位置 | 文献识别信息 |','|---|---|---|']
    for n,places in sorted(citation_locations(text).items()):
        quoted = re.search('“([^”]+)”', refs[n])
        linked = re.search(r'\[([^\]]+)\]\(https://', refs[n])
        title = quoted.group(1) if quoted else (linked.group(1) if linked and not linked.group(1).startswith('arXiv:') else refs[n].split(', ', 1)[1].rsplit(', [arXiv:', 1)[0])
        lines.append(f'| [{n}] | {"；".join(places)} | {title} |')
    return '\n'.join(lines)+'\n\n'

def check_all(here):
    here=Path(here)
    text=(here/'MANUSCRIPT.md').read_text(encoding='utf-8')
    sec=sections(text)
    manifest=json.loads((here/'translation_sync.json').read_text(encoding='utf-8'))
    assert manifest['scope']=='Reviewed version fingerprints plus structural checks; not automatic semantic equivalence.'
    assert list(manifest['sections']) == [key for key,_ in PAIRS]
    hashes = lambda s: hashlib.sha256(s.encode('utf-8')).hexdigest()
    checked_links=0
    for key,name in PAIRS:
        zh=(here/name).read_text(encoding='utf-8')
        zh_body=zh.split('\n## ',1)[1]
        assert_text_pair(sec[key],zh_body)
        entry=manifest['sections'][key]
        assert entry['chinese_file']==name
        assert entry['english_sha256']==hashes(sec[key]), f'{key}: English version changed since review'
        assert entry['chinese_sha256']==hashes(zh), f'{name}: translation changed since review'
        assert not re.search('去极化|重量|非一致|多重线性',zh_body),f'{name}: noncanonical terminology'
        for link in re.findall(r'\]\(([^)]+)\)',zh):
            if not link.startswith(('https://','http://')):
                assert (here/link).is_file(),(name,link)
                checked_links+=1
        # Figure/table captions must retain their IDs in each translation.
        en_ids=re.findall(r'^\*\*(?:Figure|Table) ([SC]?\d+)\.',sec[key],re.M)
        zh_ids=re.findall(r'^\*\*(?:图|表)\s*([SC]?\d+)[：.:]',zh_body,re.M)
        assert en_ids==zh_ids,(name,'Caption IDs',en_ids,zh_ids)
    registry=(here/'REFERENCES_MANAGEMENT_ZH.md').read_text(encoding='utf-8')
    expected=reference_index(text)
    actual=registry[registry.index('## 当前引用位置索引'):registry.index('## 正式引用清单')]
    assert actual==expected,'Reference location index stale'
    entries=dict((int(n),v) for n,v in re.findall(r'^### \[(\d+)\][^\n]*\n(.*?)(?=^### \[|\Z)',registry,re.M|re.S))
    refs=source_refs(text)
    assert set(citation_locations(text))==set(refs)==set(entries)
    for n,v in refs.items():
        assert re.search(r'\*\*文献：\*\* (.+)',entries[n]).group(1).strip()==v, f'Reference {n} differs'
    # Recompute abstract upper bound from full-grid source, not a stored expected ratio.
    with (here.parent/'results/paper1_grid/grid_ler.csv').open(encoding='utf-8-sig',newline='') as f:
        grid=list(csv.DictReader(f))
    by={(x['d'],x['p'],x['method']):x for x in grid}
    ratios=[float(x['LER'])/float(by[x['d'],x['p'],x['reference']]['LER']) for x in grid if x['is_main']=='True' and int(x['d'])<=9]
    abstract=sec['Abstract']
    bound=float(re.search(r'at most ([0-9.]+) times',abstract).group(1))
    assert max(ratios)<=bound,'Abstract upper bound violated'
    assert str(bound) in (here/'ABSTRACT_ZH.md').read_text(encoding='utf-8')
    assert 'bifurcation sweeps' not in abstract
    # Archived Appendix C remains auditable but is outside the submission.
    with (here.parent/'results/circuit_locality/summary/locality_summary.csv').open(encoding='utf-8-sig',newline='') as f:
        stats=list(csv.DictReader(f))
    for table in [(here/'archive/APPENDIX_C_ARCHIVED.md').read_text(encoding='utf-8'),
                  (here/'APPENDIX_C_ZH.md').read_text(encoding='utf-8')]:
        rows=[l for l in table.splitlines() if l.startswith(('| Surface','| 表面码','| BB '))]
        assert len(rows)==len(stats)==7
        for line,x in zip(rows,stats):
            cells=[c.strip() for c in line.split('|')[1:-1]]
            for i,field in [(1,'n_mech'),(2,'dim_spins'),(4,'n_global'),(6,'order_max')]:
                assert int(cells[i].replace(',',''))==int(x[field]),('Table C1',field)
            p50,p99=map(float,cells[3].split('/'))
            assert abs(p50-float(x['self_w_p50']))<1e-10
            assert abs(p99-float(x['self_w_p99']))<1e-10
            assert cells[5]=='—' if int(x['n_global'])==0 else int(cells[5])==int(x['global_w_max'])
    return {'formal_translation_sections':len(PAIRS),'translation_links_valid':checked_links,
            'reference_location_index_current':True,'reference_entries_identical':True,
            'display_math_images_citations_captions_synced':True,'reviewed_version_hashes_current':True,
            'small_code_max_ratio':max(ratios),'abstract_ratio_bound':bound,
            'table_C1_rows_checked_in_each_language':len(stats),
            'limitation':'Fingerprints detect later edits; structural checks do not prove semantic translation equivalence.'}
