# -*- coding: utf-8 -*-
"""从朱宁 rebuild 模型的 WCM 独立输入生成 test/A、test/B 两套待测包。

源文件只读；唯一允许的内容改动是：WCM Bin 本构卡换成 UMAT 材料卡与
精简 SDV 输出、UVARM 场请求换成编号 SDV、场输出为 LE/PE/PEEQ/PEMAG/S
与 U、*Static 增量、B 版复材截面的 Enhanced 沙漏控制。网格、铺层、
约束、载荷、取向一律不动，末尾用保护摘要证明。

A 版按用户要求不加 *Section Controls；B 版（能量线性软化）加
hourglass=ENHANCED。两套材料卡与各自 PROPS 项数（A 25 / B 30）配套。

用法：python tools/prepare_test_cases.py [--check] [--overwrite]
--check 只重建并比对字节，不写任何文件。
"""
import argparse
import csv
import hashlib
import json
import math
import re
from collections import Counter
from pathlib import Path
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
from convert_wcm_inp import blocks, data_end, material_groups, transform
from model import output_indices, sdv_labels, state_count, validate_props

ROOT = Path(__file__).resolve().parents[1]          # wcm_hashin_umat/
REPO = ROOT.parent
CONTROL = 'HG_UMAT'
# 朱宁 rebuild 基准，旧 10deg 的 Job-2 已废弃，不再作为任何输入。
SOURCE_INP = REPO/'geometry_zhuning_joint_rebuild'/'WCM_Job.inp'
SOURCE_ANG = REPO/'geometry_zhuning_joint_rebuild'/'WCM_Zhuning_Joint_Rebuilt_WindAngles.ang'
TANK_PART = 'tank-1'
STATIC = [0.01, 1.0, 1e-20, 0.1]
FIELD_ELEMENT = ['LE', 'PE', 'PEEQ', 'PEMAG', 'S', 'SDV']
EXPECTED = {'A': {'elements': 51558, 'bins': 135, 'controls': False},
            'B': {'elements': 51558, 'bins': 135, 'controls': True}}
CASES = {'A': {'dir': REPO/'test'/'A_stiffness_reduction', 'inp': 'Rebuild_PDM_A.inp',
               'umat': 'hashin_constant_wcm.for'},
         'B': {'dir': REPO/'test'/'B_energy_softening', 'inp': 'Rebuild_PDM_B.inp',
               'umat': 'hashin_energy_wcm.for'}}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def content(block):
    return block.text[:data_end(block) - block.start].rstrip()


def patch(text, edits):
    edits = sorted(edits)
    require(all(a[1] <= b[0] for a, b in zip(edits, edits[1:])), '替换区间重叠')
    pieces, position = [], 0
    for start, end, replacement in edits:
        pieces.extend((text[position:start], replacement))
        position = end
    pieces.append(text[position:])
    return ''.join(pieces)


def adapt_requests(text):
    """UVARM 场请求换成编号 SDV，输出量与步长按用户口径固定。"""
    items = blocks(text)
    fields = [b for b in items if b.key == 'element output']
    nodes = [b for b in items if b.key == 'node output']
    statics = [b for b in items if b.key == 'static']
    require(len(fields) == len(nodes) == len(statics) == 1,
            '仅接受单步、单场输出块的 WCM 独立导出；当前结构不符')
    require(re.search(r'\bUVARM\d*\b', fields[0].text, re.I), '源场请求里没有 UVARM')
    increments = numbers(statics[0])
    require(increments == STATIC, f'*Static 增量与要求不一致：{increments}')
    edits = [(fields[0].start, data_end(fields[0]),
              '*Element Output, directions=YES\n' + ', '.join(FIELD_ELEMENT) + '\n'),
             (nodes[0].start, data_end(nodes[0]), '*Node Output\nU\n')]
    for b in items:
        head = b.text.splitlines()[0].lower()
        if b.key == 'output' and 'history' in head:
            edits.append((b.start, data_end(b), ''))
        if b.key in ('contact output', 'energy output'):
            edits.append((b.start, data_end(b), ''))
    return patch(text, edits)


def data_lines(block):
    return [s.strip() for s in block.text.splitlines()[1:]
            if s.strip() and not s.lstrip().startswith('**')]


def numbers(block):
    tokens = [x.strip() for line in data_lines(block) for x in line.split(',') if x.strip()]
    return [float(t.replace('D', 'E').replace('d', 'e')) for t in tokens]


def add_controls(text, targets):
    """仅 B 版：给复材截面挂 Enhanced 沙漏控制。"""
    items = blocks(text)
    require(not any(b.key in ('section controls', 'hourglass stiffness') for b in items),
            '源文件已有沙漏设置，拒绝叠加')
    edits = []
    for b in items:
        if b.key == 'solid section' and b.options.get('material', '').lower() in targets:
            require('controls' not in b.options, '目标截面已有控制')
            line = b.text.splitlines()[0]
            edits.append((b.start, b.start + len(line), line + ', controls=' + CONTROL))
    require(edits, '未找到复材截面')
    first_part = next(b for b in items if b.key == 'part')
    edits.append((first_part.start, first_part.start,
                  '*Section Controls, name=' + CONTROL + ', hourglass=ENHANCED\n'))
    return patch(text, edits), len(edits) - 1


def tank_topology(text):
    """解析复材 part 的单元/集合/截面，用于逐单元映射审计。"""
    items = blocks(text)
    active = False
    elements, sets, sections, orientations, distributions = {}, {}, [], {}, {}
    for b in items:
        if b.key == 'part':
            active = b.options.get('name', '').lower() == TANK_PART
        elif b.key == 'end part':
            active = False
        if not active:
            continue
        lines = data_lines(b)
        if b.key == 'element':
            kind = b.options['type'].upper()
            require(kind in ('C3D8', 'C3D8R', 'C3D6'), '不支持的复材单元：' + kind)
            for line in lines:
                vals = [int(x) for x in line.split(',') if x.strip()]
                require(len(vals) == (7 if kind == 'C3D6' else 9), '连接行格式不符')
                require(vals[0] not in elements, '单元重复')
                elements[vals[0]] = kind
                if 'elset' in b.options:
                    sets.setdefault(b.options['elset'].lower(), []).append(str(vals[0]))
        elif b.key == 'elset':
            name = b.options['elset'].lower()
            values = sets.setdefault(name, [])
            for line in lines:
                tokens = [x.strip() for x in line.split(',') if x.strip()]
                if 'generate' in b.options:
                    nums = [int(x) for x in tokens]
                    require(len(nums) == 3 and nums[2] > 0, '无效 generate')
                    values.extend(str(x) for x in range(nums[0], nums[1] + 1, nums[2]))
                else:
                    values.extend(tokens)
        elif b.key == 'solid section':
            sections.append(b)
        elif b.key == 'orientation':
            orientations[b.options['name'].lower()] = lines
        elif b.key == 'distribution':
            require(b.options.get('location', '').lower() == 'element', '分布不是逐单元')
            values = {}
            for line in lines:
                label, angle = [x.strip() for x in line.split(',')]
                values[int(label) if label else None] = float(angle)
            distributions[b.options['name'].lower()] = values
    require(elements and sections, '复材 part 单元或截面缺失')

    cache = {}

    def resolve(name, trail=()):
        name = name.lower()
        require(name not in trail and name in sets, '集合引用缺失/循环：' + name)
        if name not in cache:
            result = set()
            for token in sets[name]:
                result.add(int(token)) if token.isdigit() else result.update(
                    resolve(token, trail + (name,)))
            cache[name] = result
        return cache[name]

    return elements, sections, orientations, distributions, resolve


def read_angles(path, labels):
    rows = {}
    for row in csv.reader(path.read_text(encoding='ascii').splitlines()):
        if not row:
            continue
        require(len(row) == 7, '角度文件应为七列')
        require(tuple(int(x) for x in row[1:6]) == (1, 1, 1, 1, 1), '角度索引与本基准不符')
        label = int(row[0])
        angle = float(row[6])
        require(label not in rows and math.isfinite(angle), '角度重复/非有限')
        require(-1e-12 <= angle <= math.pi/2 + 1e-12, '弧度角超范围')
        rows[label] = angle
    require(set(rows) == set(labels), 'ang 与当前复材网格不是逐单元一一对应')
    return rows


def audit(deck, records, variant, want_controls):
    """逐单元截面/角度一致性；返回审计摘要。"""
    materials = {r['material'].lower(): r for r in records}
    theta_index = 22 if variant == 'A' else 25
    elements, sections, orientations, distributions, resolve = tank_topology(deck)
    assignment = {}
    for section in sections:
        material = section.options.get('material', '').lower()
        require(material in materials, '复材 part 里存在非目标材料')
        require(section.options.get('controls') == (CONTROL if want_controls else None),
                '复材截面的沙漏控制引用与本版要求不符')
        require(section.options.get('orientation'), '复材截面缺取向')
        for label in resolve(section.options['elset']):
            require(label in elements and label not in assignment,
                    '截面漏引用/重复赋值：' + str(label))
            assignment[label] = (material, section.options['orientation'])
    require(set(assignment) == set(elements), '复材单元截面覆盖不全')
    ang = read_angles(SOURCE_ANG, elements)
    deltas = []
    addrot = {}
    for label in sorted(elements):
        material, orientation = assignment[label]
        lines = orientations[orientation.lower()]
        require(len(lines) == 2, '取向数据行数异常')
        axis, dist_name = [x.strip() for x in lines[1].split(',')]
        require(axis == '2', 'WCM AddRot 旋转轴改变')
        table = distributions[dist_name.lower()]
        addrot[label] = table.get(label, table.get(None))
        require(addrot[label] is not None, 'AddRot 未覆盖单元')
        record = materials[material]
        require(abs(record['props'][theta_index] - record['angle']) <= 1e-10,
                'Beta 与 PROPS 角度不一致')
        deltas.append(record['angle'] - math.degrees(ang[label]))
    return {'composite_elements': len(elements), 'material_bins': len(materials),
            'sections_checked': len(sections),
            'element_types': dict(Counter(elements.values())),
            'bin_minus_ang_deg_min': min(deltas), 'bin_minus_ang_deg_max': max(deltas),
            'addrot_defined': len(addrot)}


def protected_digest(text, targets):
    """保护摘要：许可改动（目标 Bin 材料卡、输出/步长请求、B 的沙漏控制）
    之外的字节必须逐字一致，用来证明网格/铺层/约束/载荷未被触碰。"""
    items = blocks(text)
    excluded = set()
    for key, (_, heading, group) in material_groups(items).items():
        if key in targets:
            excluded.update(b.start for b in [heading] + group)
    rows = []
    for b in items:
        if b.start in excluded or b.key in ('static', 'node output', 'element output',
                                           'contact output', 'energy output', 'output'):
            continue
        if b.key == 'section controls' and b.options.get('name') == CONTROL:
            continue
        value = content(b)
        # 注释不影响模型身份；关键字与数据的其余字节保留。
        value = '\n'.join(s for s in value.splitlines() if not s.lstrip().startswith('**'))
        if b.key == 'solid section':
            value = re.sub(r', controls=' + CONTROL + r'(?=\n|$)', '', value)
        rows.append(value)
    return digest('\n'.join(rows).encode('latin-1'))


def verify_deck(deck, variant, records):
    items = blocks(deck)
    groups = material_groups(items)
    count = state_count(variant)
    for record in records:
        group = groups[record['material'].lower()][2]
        user = [b for b in group if b.key == 'user material']
        depvar = [b for b in group if b.key == 'depvar']
        require(len(user) == len(depvar) == 1 and 'unsymm' in user[0].options, '本构卡结构不符')
        props = user[0].numbers()
        validate_props(variant, props)
        require(int(user[0].options['constants']) == len(props), '常数声明与实参不符')
        require(len(props) == (25 if variant == 'A' else 30), 'PROPS 项数错误')
        require(np.allclose(props, record['props'], rtol=1e-11, atol=1e-12), '写出参数精度异常')
        require(data_lines(depvar[0]) == [f'{count},'], '状态数不符或仍带自定义 SDV 标签')
        require(not any(b.key in ('elastic', 'fail stress', 'fail strain',
                                  'user output variables') for b in group),
                '目标材料残留 WCM 本构/输出')
    require(not any(b.key == 'user output variables' for b in items), '残留 UVARM 定义')
    fields = [content(b) for b in items if b.key == 'element output']
    require(not re.search(r'\bUVARM\d*\b', '\n'.join(fields), re.I), '残留 UVARM 场请求')
    body = '\n'.join(fields).upper()
    want = [f'SDV{i}' for i in output_indices(variant)]
    for name in FIELD_ELEMENT[:-1] + ['S'] + want:
        require(re.search(r'\b' + name + r'\b', body), '场输出缺少 ' + name)
    for name in re.findall(r'\bSDV(\d+)\b', body):
        require(int(name) <= len(want), '请求了未导出的 SDV' + name)
    nodes = [content(b) for b in items if b.key == 'node output']
    require(all(b.splitlines()[-1].strip().upper() == 'U' for b in nodes), '节点输出应为 U')
    statics = numbers(next(b for b in items if b.key == 'static'))
    require(statics == STATIC, '步长与要求不符')
    return {'state_variables': count, 'output_sdvs': want, 'static': statics,
            'element_output': FIELD_ELEMENT[:-1] + ['S'] + want, 'node_output': ['U'],
            'material_cards': len(records),
            'props_constants': 25 if variant == 'A' else 30}


def build(variant, config, source):
    case = CASES[variant]
    adapted = adapt_requests(source)
    deck, records = transform(adapted, variant, config)
    targets = {r['material'].lower() for r in records}
    want_controls = EXPECTED[variant]['controls']
    controls = 0
    if want_controls:
        deck, controls = add_controls(deck, targets)
    checks = verify_deck(deck, variant, records)
    audit_info = audit(deck, records, variant, want_controls)
    require(audit_info['composite_elements'] == EXPECTED[variant]['elements'], '单元数与基准不符')
    require(audit_info['material_bins'] == EXPECTED[variant]['bins'], 'Bin 数与基准不符')
    require(digest(pack_umat(variant)) == digest(
        (ROOT/'dist'/case['umat']).read_bytes()), '打包源码与 dist 不一致')
    info = {**checks, **audit_info, 'hourglass_enhanced_sections': controls,
            'max_relative_elastic_error': max(r['relative_elastic_error'] for r in records),
            'sdv_labels': dict(enumerate(sdv_labels(variant), 1))}
    return deck, info, records


def pack_umat(variant):
    """模板 + 核心拼装，并与 dist/ 的构建产物保持逐字节一致。"""
    template = (ROOT/'src'/'umat_entry.for.in').read_text(encoding='utf-8')
    core = (ROOT/'src'/'wcm_core.for').read_text(encoding='utf-8')
    mode = 1 if variant == 'A' else 2
    name = 'A（直接刚度折减）' if variant == 'A' else 'B（能量线性软化）'
    from build import props_layout
    text = (template.replace('@MODE@', str(mode))
            .replace('@VERSION@', name)
            .replace('C @PROPS_LAYOUT@', props_layout(mode)) + '\n' + core)
    for number, line in enumerate(text.splitlines(), 1):
        if not line or line[0] in 'Cc*!':
            continue
        require('\t' not in line and len(line) <= 72 and line.isascii(),
                'Fortran 固定格式不合法：' + str(number))
    hashin = core.split('SUBROUTINE WCM_HASHIN', 1)[1].split('SUBROUTINE WCM_EQ', 1)[0]
    require(hashin.count('(S(5)/P(14))**2') == 3 and '(S(5)/P(15))' not in hashin,
            '朱宁 tau13/S12 判据不符合当前确认版本')
    require('@MODE@' not in text and '@VERSION@' not in text, '构建占位符残留')
    return text.encode('utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='只重建比对，不写文件')
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--config', type=Path, default=ROOT/'config'/'project.json')
    args = parser.parse_args()
    raw = SOURCE_INP.read_bytes()
    source = raw.decode('latin-1').replace('\r\n', '\n')
    config = json.loads(args.config.read_text(encoding='utf-8'))
    require(config.get('weight_plus') == 0.5, '本轮仅支持 WCM 平衡双族')
    require(config.get('update') == 0, '更新算法变更须另立对照，不能静默修改')
    require((config.get('fiber_cap'), config.get('matrix_cap')) == (1.0, 1.0),
            'B 刚度折减上限应为 1，损伤允许演化到 1')
    report = {'source_inp': str(SOURCE_INP), 'source_sha256': digest(raw),
              'source_ang_sha256': digest(SOURCE_ANG.read_bytes()),
              'line_endings': 'CRLF' if b'\r\n' in raw else 'LF',
              'baseline': 'geometry_zhuning_joint_rebuild（旧 10deg Job-2 已废弃）',
              'cases': {}}
    for variant in ('A', 'B'):
        deck, info, records = build(variant, config, source)
        text = deck.replace('\n', '\r\n') if b'\r\n' in raw else deck
        case = CASES[variant]
        artifacts = {case['dir']/case['inp']: text.encode('latin-1'),
                     case['dir']/case['umat']: pack_umat(variant)}
        for path, data in artifacts.items():
            require(path.parent.is_dir(), '目标套件目录缺失：' + str(path.parent))
            if args.check:
                require(path.is_file() and path.read_bytes() == data,
                        '产物与当前生成规则不同：' + str(path))
            else:
                if path.exists() and not args.overwrite:
                    raise FileExistsError('已存在，需 --overwrite：' + str(path))
                path.write_bytes(data)
        report['cases'][variant] = {
            'inp': str(case['dir']/case['inp']), 'umat': case['umat'],
            'inp_sha256': digest(artifacts[case['dir']/case['inp']]),
            'umat_sha256': digest(artifacts[case['dir']/case['umat']]),
            'protected_source_sha256': protected_digest(source,
                                                        {r['material'].lower() for r in records}),
            'protected_output_sha256': protected_digest(deck,
                                                        {r['material'].lower() for r in records}),
            'materials': [{'material': r['material'], 'angle_deg': r['angle'],
                           'ply': r['ply'], 'elastic_error': r['relative_elastic_error']}
                          for r in records],
            **{k: v for k, v in info.items() if k != 'sdv_labels'}}
        require(report['cases'][variant]['protected_source_sha256'] ==
                report['cases'][variant]['protected_output_sha256'],
                '许可范围外的模型内容发生变化')
    # 两套之间除 Bin 材料卡与 B 的沙漏控制外必须逐字一致。
    require(report['cases']['A']['protected_output_sha256'] ==
            report['cases']['B']['protected_output_sha256'], 'A/B 非材料内容出现差异')
    report['sdv_labels'] = {v: dict(enumerate(sdv_labels(v), 1)) for v in ('A', 'B')}
    report['written'] = not args.check
    report['note'] = ('A 版按用户口径不加沙漏控制；PE/PEEQ/PEMAG 由 UMAT 不返回非弹性应变，'
                      '预期为零值，仅供与 WCM 弹性基准对照。')
    out = ROOT/'validation'/'test_cases.json'
    if args.check:
        previous = json.loads(out.read_text(encoding='utf-8'))
        for key in ('source_sha256', 'cases'):
            require(previous[key] == report[key], '校验结果与已交付产物不一致：' + key)
    else:
        out.write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    summary = {v: {k: report['cases'][v][k] for k in
                   ('inp_sha256', 'umat_sha256', 'material_cards', 'props_constants',
                    'state_variables', 'output_sdvs', 'hourglass_enhanced_sections',
                    'composite_elements', 'bin_minus_ang_deg_min', 'bin_minus_ang_deg_max',
                    'max_relative_elastic_error')} for v in ('A', 'B')}
    print(json.dumps({'action': 'check' if args.check else 'prepare',
                      'report': str(out), **summary}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
