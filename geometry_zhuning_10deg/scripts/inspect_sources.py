# -*- coding: utf-8 -*-
"""只读核对朱宁论文与旧 INP；证据输出仅写入本目录。"""
from pathlib import Path
import hashlib
import json
import math
import fitz

ROOT = Path(__file__).resolve().parent.parent  # 脚本位于 scripts\，数据根目录为上级
PROJECT = ROOT.parent
OUT = ROOT / 'evidence'
OUT.mkdir(exist_ok=True)
pdf = PROJECT / 'references' / 'REF3_朱宁-2023_浙江大学硕士论文_模型铺层与材料参数来源.pdf'
with fitz.open(str(pdf)) as doc:
    texts = []
    for number in (31, 32, 33, 46, 47, 48, 52, 57):
        page = doc[number - 1]
        texts.append('===== PDF PAGE {} =====\n{}'.format(number, page.get_text()))
        if number in (33, 48, 52):
            page.get_pixmap(matrix=fitz.Matrix(2.5, 2.5)).save(str(OUT / ('zhuning_page_{}.png'.format(number))))
    (OUT / 'zhuning_geometry_extract.txt').write_text('\n\n'.join(texts), encoding='utf-8')

inp = PROJECT / 'Job-hashinnewcdm.inp'
parts = {}
part = None
nodes = False
with inp.open(encoding='utf-8', errors='replace') as stream:
    for raw in stream:
        line = raw.strip()
        if not line or line.startswith('**'):
            continue
        if line.startswith('*'):
            nodes = False
            if line.lower().startswith('*part,'):
                name = line.split('name=', 1)[1].strip()
                part = {'node_count': 0, 'min_xyz': [float('inf')] * 3,
                        'max_xyz': [-float('inf')] * 3, 'radius_min': float('inf'),
                        'radius_max': 0.0, 'angle_min': float('inf'),
                        'angle_max': -float('inf'), 'sections': []}
                parts[name] = part
            elif line.lower().startswith('*end part'):
                part = None
            elif line.lower() == '*node' and part is not None:
                nodes = True
            elif line.lower().startswith('*solid section') and part is not None:
                if 'WCM_' not in line:
                    part['sections'].append(line)
            continue
        if nodes:
            values = [float(x) for x in line.split(',')[1:4]]
            x, y, z = values
            part['node_count'] += 1
            part['min_xyz'] = [min(a, b) for a, b in zip(part['min_xyz'], values)]
            part['max_xyz'] = [max(a, b) for a, b in zip(part['max_xyz'], values)]
            radius = math.hypot(x, z)
            angle = math.degrees(math.atan2(z, x))
            part['radius_min'] = min(part['radius_min'], radius)
            part['radius_max'] = max(part['radius_max'], radius)
            part['angle_min'] = min(part['angle_min'], angle)
            part['angle_max'] = max(part['angle_max'], angle)

report = {'paper': str(pdf), 'paper_sha256': hashlib.sha256(pdf.read_bytes()).hexdigest(),
          'inp': str(inp), 'inp_sha256': hashlib.sha256(inp.read_bytes()).hexdigest(),
          'parts': parts,
          'old_cad_candidates': [r'G:\cylinder\qiping\model-sw\zhuning-cylinder.STEP',
                                 r'G:\cylinder\qiping\model-cad\cylinder2-zhuning.dxf',
                                 r'G:\cylinder\qiping\model-cad\zhuning-cylinder.DXF'],
          'note': '历史 CAD 仅为辅助依据；不得替代论文已明确给出的尺寸。'}
(OUT / 'source_audit.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(parts, ensure_ascii=False, indent=2))
print('EVIDENCE_READY', str(OUT))

# 提取原图局部用于几何比例核对，不将截图量测值当作论文标注值。
with fitz.open(str(pdf)) as doc:
    doc[32].get_pixmap(matrix=fitz.Matrix(8, 8), clip=fitz.Rect(418, 630, 460, 701)).save(str(OUT / 'zhuning_boss_detail.png'))

# 仅解析 DXF 的文本实体参数，不修改历史 CAD。
dxf = Path(r'G:\cylinder\qiping\model-cad\cylinder2-zhuning.dxf')
rows = dxf.read_text(encoding='gb18030', errors='replace').splitlines()
pairs = [(rows[i].strip(), rows[i + 1].strip()) for i in range(0, len(rows) - 1, 2)]
entities, current, in_entities = [], None, False
for code, value in pairs:
    if code == '2' and value == 'ENTITIES':
        in_entities = True
    if not in_entities:
        continue
    if code == '0':
        if current is not None:
            entities.append(current)
        if value == 'ENDSEC':
            break
        current = {'type': value, 'groups': {}}
    elif current is not None:
        current['groups'].setdefault(code, []).append(value)
selected = [e for e in entities if e['groups'].get('8') == ['0']]
(OUT / 'legacy_dxf_entities.json').write_text(json.dumps(selected, ensure_ascii=False, indent=2), encoding='utf-8')
for e in selected:
    g = e['groups']
    if e['type'] in ('LINE', 'ARC', 'ELLIPSE'):
        print(e['type'], {k: g[k] for k in ('5', '10', '20', '11', '21', '40', '41', '42', '50', '51', '230') if k in g})
