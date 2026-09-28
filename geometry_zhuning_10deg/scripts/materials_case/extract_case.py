# -*- coding: utf-8 -*-
"""只读提取原案例三张材料卡及筒身铺层，不修改任何历史输入文件。"""
from pathlib import Path
import hashlib
import json
import math
import re

ROOT = Path(__file__).resolve().parent.parent.parent  # 数据根目录：脚本位于 scripts\materials_case\
PROJECT = ROOT.parent
NAMES = ('al-6061', 't-700', 'Liner')
SOURCES = [PROJECT / 'Job-hashinnewcdm.inp', PROJECT / 'Job-hashinnewconstant.inp']


def keywords(text):
    result = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith('**'):
            continue
        if line.startswith('*'):
            tokens = [q.strip() for q in line[1:].split(',')]
            opts = {}
            for token in tokens[1:]:
                key, sep, value = token.partition('=')
                opts[key.lower()] = value if sep else True
            result.append({'keyword': tokens[0].lower(), 'options': opts,
                           'header': line, 'rows': []})
        else:
            assert result, '数据行之前缺少关键字'
            result[-1]['rows'].append(line)
    return result


def material_cards(blocks):
    result, active = {}, None
    supported = {'density', 'elastic', 'plastic', 'fail strain', 'fail stress'}
    for block in blocks:
        if block['keyword'] == 'material':
            active = block['options']['name']
            if active in NAMES:
                result[active] = {'header': block['header'], 'properties': []}
        elif active in NAMES:
            if block['keyword'] not in supported:
                active = None
                continue
            prop = dict(block)
            prop['values'] = [[float(x.strip()) for x in row.split(',') if x.strip()]
                              for row in block['rows']]
            result[active]['properties'].append(prop)
    assert set(result) == set(NAMES)
    return result


def labels(block):
    values = [int(q.strip()) for row in block['rows'] for q in row.split(',') if q.strip()]
    if 'generate' in block['options']:
        assert len(values) % 3 == 0
        return set(v for i in range(0, len(values), 3)
                   for v in range(values[i], values[i + 1] + 1, values[i + 2]))
    return set(values)


texts = [p.read_text(encoding='utf-8') for p in SOURCES]
blocks = keywords(texts[0])
cards = material_cards(blocks)
assert cards == material_cards(keywords(texts[1])), '两个历史案例的基础材料卡不一致'

nodes, elements, nsets, elsets, sections = {}, {}, {}, {}, []
in_tank = False
for block in blocks:
    key, opts = block['keyword'], block['options']
    if key == 'part':
        in_tank = opts['name'] == 'Tank-1'
    elif key == 'end part':
        in_tank = False
    elif in_tank:
        if key == 'node':
            for row in block['rows']:
                data = [q.strip() for q in row.split(',') if q.strip()]
                nodes[int(data[0])] = tuple(float(q) for q in data[1:4])
        elif key == 'element':
            for row in block['rows']:
                data = [int(q.strip()) for q in row.split(',') if q.strip()]
                elements[data[0]] = data[1:]
        elif key == 'nset' and re.fullmatch(r'Layer\d+', opts['nset']):
            nsets.setdefault(opts['nset'], set()).update(labels(block))
        elif key == 'elset':
            elsets.setdefault(opts['elset'], set()).update(labels(block))
        elif key == 'solid section':
            sections.append((opts['elset'], opts['material']))

beta_by_material = {}
last_beta = None
for line in texts[0].splitlines():
    match = re.search(r'\*\*.*Beta\s*=\s*([\d.]+)', line)
    if match:
        last_beta = float(match[1])
    elif line.lower().startswith('*material,') and 'WCM_' in line:
        beta_by_material[line.split('name=', 1)[1].strip()] = last_beta

center_elements = {label for label, connectivity in elements.items()
                   if any(abs(nodes[n][1]) < 1e-8 for n in connectivity)
                   and max(nodes[n][1] for n in connectivity) > 1e-8
                   and min(nodes[n][1] for n in connectivity) >= -1e-8}
materials_at_center = {}
for set_name, material in sections:
    for label in elsets[set_name] & center_elements:
        if label in materials_at_center:
            assert materials_at_center[label] == material
        materials_at_center[label] = material

layup = []
for name in sorted(nsets):
    center_nodes = [nodes[n] for n in nsets[name] if abs(nodes[n][1]) < 1e-8]
    radii = [math.hypot(q[0], q[2]) for q in center_nodes]
    angles = {beta_by_material[materials_at_center[e]] for e in elsets[name] & center_elements}
    assert len(angles) == 1, (name, angles)
    lo, hi = min(radii), max(radii)
    thickness = round(hi - lo, 3)
    count = round(thickness / .3)
    assert abs(thickness - count * .3) < .002
    layup.append({'group': name, 'angle_deg': angles.pop(), 'inner_radius_mm': round(lo, 3),
                  'outer_radius_mm': round(hi, 3), 'thickness_mm': thickness,
                  'ply_count_at_0_3_mm': count})
assert len(layup) == 12
print('原案例筒身层组核对：', json.dumps(layup, ensure_ascii=False))
actual_count = sum(q['ply_count_at_0_3_mm'] for q in layup)
# 论文表3.2（PDF第52页）单独保存，不能用旧INP的74层冒充朱宁72层方案。
paper_angles = (15, 90, 15, 90, 20, 90, 25, 90, 35, 90, 40, 90)
paper_counts = (6, 6, 6, 6, 8, 8, 8, 8, 4, 4, 4, 4)
paper_layup = [{'group': 'Layer{:02d}'.format(i + 1), 'angle_deg': angle,
                'ply_count': count, 'thickness_mm': round(count * .3, 3)}
               for i, (angle, count) in enumerate(zip(paper_angles, paper_counts))]
assert all(abs(layup[i]['outer_radius_mm'] - layup[i + 1]['inner_radius_mm']) < .002
           for i in range(len(layup) - 1))

report = {'sources': [{'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in SOURCES],
          'three_cards_identical_between_cases': True, 'materials': cards,
          'layup_inner_to_outer': layup, 'nominal_ply_thickness_mm': .3,
          'actual_total_ply_count': actual_count,
          'actual_total_thickness_mm': round(sum(q['thickness_mm'] for q in layup), 3),
          'zhuning_table_3_2_inner_to_outer': paper_layup,
          'zhuning_total_ply_count': sum(paper_counts),
          'zhuning_total_thickness_mm': round(sum(paper_counts) * .3, 3),
          'thickness_basis': '原INP的Tank-1中Layer01至Layer12在Y=0的节点半径差，并按对应单元材料前的Beta注释核对角度。层数按朱宁表3.2的0.3 mm单层计。',
          'half_cylinder_length_mm': 641.0,
          'notes': ['Liner原卡未定义Density，按原案例保持未定义。',
                    't-700为原始单层正交各向异性卡；不复制WCM等效Bin，也不预设尚在修改的UMAT参数。',
                    'Fail Stress和Fail Strain按旧卡复制，不代表已经配置Hashin渐进损伤演化。',
                    '螺旋组的层数为正负角合计，例如±15°的6层厚度为1.8 mm，不是3.6 mm。']}
ROOT.mkdir(exist_ok=True)
(ROOT / 'case_materials_and_layup.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
lines = ['*Heading', '** 仅包含旧案例材料卡，不包含网格、截面指派或分析步。']
for name in NAMES:
    card = cards[name]
    lines.append(card['header'])
    for prop in card['properties']:
        lines.extend([prop['header']] + prop['rows'])
(ROOT / 'case_material_cards.inp').write_text('\n'.join(lines) + '\n', encoding='utf-8')
print(json.dumps({'layup': layup, 'material_names': list(cards), 'checks': report['notes']}, ensure_ascii=False, indent=2))
