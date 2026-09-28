# -*- coding: utf-8 -*-
"""保留WCM模型，只替换指定Bin材料。源文件只读，输出限于新仓库。

支持当前独立inp和WCM直接导出的独立inp。为避免漏转嵌套材料，
遇到*Include主动拒绝，请先导出独立inp；不猜测不完整的材料映射。
既可读取旧A/B用户材料，也可读取原始WCM工程常数。原9项用于
验证等效弹性兼容性，输出前9项则来自原始单层材料t-700。
"""
from pathlib import Path
from dataclasses import dataclass
import argparse
import csv
import fnmatch
import hashlib
import json
import re
import numpy as np
from model import elastic, make_props, smeared_stiffness, state_count, output_indices

ROOT = Path(__file__).resolve().parents[1]
KEYWORD = re.compile(r'^\*(?!\*)([^\r\n]*)', re.M)
BETA = re.compile(r'^\*\*\s*\*+\s*Beta\s*=\s*([-+0-9.eEdD]+)', re.M | re.I)
MATERIAL_KEYS = {'density', 'elastic', 'plastic', 'depvar', 'user material',
                 'fail stress', 'fail strain', 'user output variables', 'expansion',
                 'damping', 'specific heat', 'conductivity'}


@dataclass
class Block:
    start: int
    end: int
    key: str
    options: dict
    text: str

    def numbers(self):
        lines = self.text.splitlines()[1:]
        tokens = [x.strip() for line in lines if not line.startswith('**')
                  for x in line.split(',') if x.strip()]
        return [float(x.replace('D', 'E').replace('d', 'e')) for x in tokens]


def blocks(text):
    hits = list(KEYWORD.finditer(text))
    result = []
    for i, hit in enumerate(hits):
        tokens = [s.strip() for s in hit[1].split(',')]
        opts = {}
        for token in tokens[1:]:
            pair = token.split('=', 1)
            opts[pair[0].lower()] = pair[1].strip() if len(pair) == 2 else True
        end = hits[i+1].start() if i+1 < len(hits) else len(text)
        result.append(Block(hit.start(), end, tokens[0].lower(), opts, text[hit.start():end]))
    return result


def material_groups(items):
    groups = {}
    for i, block in enumerate(items):
        if block.key != 'material':
            continue
        name = block.options.get('name')
        if not name or name.lower() in groups:
            raise ValueError('重复或无名称材料：' + str(name))
        group = []
        for sub in items[i+1:]:
            if sub.key not in MATERIAL_KEYS:
                break
            group.append(sub)
        groups[name.lower()] = (name, block, group)
    return groups


def read_stiffness(group):
    candidates = [b for b in group if b.key in ('elastic', 'user material')]
    if len(candidates) != 1:
        raise ValueError('材料必须只有一项弹性或用户本构定义')
    b = candidates[0]
    if b.key == 'elastic' and b.options.get('type', '').lower() != 'engineering constants':
        raise ValueError('仅支持无温度/场依赖的工程常数')
    values = b.numbers()
    if 'dependencies' in b.options or (b.key == 'elastic' and len(values) != 9):
        raise ValueError('不支持场/温度依赖材料')
    if b.key == 'user material':
        if len(values) not in (22, 25) or int(b.options['constants']) != len(values):
            raise ValueError('仅接受旧22/25项用户材料，不可重复转换新卡')
    return values, b


def data_end(block):
    """仅替换关键词与数据行，保留尾随Beta等注释，防止误吞下一材料角度。"""
    lines = block.text.splitlines(keepends=True)
    while lines and (not lines[-1].strip() or lines[-1].lstrip().startswith('**')):
        lines.pop()
    return block.start + sum(len(s) for s in lines)


def select_sdv_outputs(text, variant):
    """只请求编号汇总量，其他场变量和尾随注释保持不变。"""
    edits = []
    for block in blocks(text):
        if block.key != 'element output':
            continue
        tokens = [v.strip() for line in block.text.splitlines()[1:]
                  if not line.lstrip().startswith('**')
                  for v in line.split(',') if v.strip()]
        is_sdv = lambda value: re.fullmatch(r'SDV(?:\d+|_\w+)?', value, re.I)
        if not any(is_sdv(v) for v in tokens):
            continue
        other = [v for v in tokens if not is_sdv(v)]
        lines = [block.text.splitlines()[0]]
        if other:
            lines.append(', '.join(other))
        names = [f'SDV{i}' for i in output_indices(variant)]
        lines += [', '.join(names[i:i+8]) for i in range(0, len(names), 8)]
        edits.append((block.start, data_end(block), '\n'.join(lines)+'\n'))
    for start, end, replacement in reversed(edits):
        text = text[:start] + replacement + text[end:]
    return text


def transform(text, variant, config):
    items = blocks(text)
    if any(b.key == 'include' for b in items):
        raise ValueError('请先导出无*Include的独立inp；拒绝可能漏转的输入')
    groups = material_groups(items)
    changes, records = [], []
    angle_hits = list(BETA.finditer(text))
    targets = set()
    for key, (name, heading, group) in groups.items():
        matching = [v for pattern, v in config['materials'].items()
                    if fnmatch.fnmatchcase(name.lower(), pattern.lower())]
        if not matching:
            continue
        if len(matching) != 1:
            raise ValueError('多重材料映射：' + name)
        mapping = matching[0]
        if any(b.key == 'depvar' and b.numbers() in ([20.0], [62.0])
               for b in group):
            raise ValueError('重复转换：已有新版 UMAT/SDV 材料卡：' + name)
        old_values, constitutive = read_stiffness(group)
        ply_key = mapping['ply_material'].lower()
        if ply_key not in groups:
            raise ValueError('缺少原始单层材料：' + ply_key)
        ply_values, ply_block = read_stiffness(groups[ply_key][2])
        if ply_block.key != 'elastic' or len(ply_values) != 9:
            raise ValueError('原始单层必须为9项工程常数')
        elastic(ply_values)
        preceding = [m for m in angle_hits if m.end() <= heading.start]
        inferred = None
        if preceding:
            hit = preceding[-1]
            gap = text[hit.end():heading.start]
            if not KEYWORD.search(gap):
                inferred = float(hit[1].replace('D', 'E').replace('d', 'e'))
        explicit = config.get('angles', {}).get(name)
        if explicit is not None and inferred is not None and abs(explicit-inferred) > 1e-10:
            raise ValueError('Beta与显式角度映射冲突：' + name)
        angle = inferred if explicit is None else float(explicit)
        if angle is None or not 0 <= angle <= 90:
            raise ValueError('缺少可靠Beta角度或超范围：' + name)
        # 用原WCM参数验证均匀化模型，而不是反推单层刚度。
        reference = elastic(old_values[:9])
        predicted = smeared_stiffness(ply_values, angle)
        error = float(np.max(np.abs(reference-predicted))/np.max(np.abs(reference)))
        if error > config.get('elastic_tolerance', 1e-4):
            raise ValueError(f'{name}的WCM弹性不兼容：{error:.6g}')
        # 旧卡的强度/损伤参数优先，保证本次不擅自改变旧参数。
        if constitutive.key == 'user material':
            expected = 22 if variant == 'A' else 25
            if len(old_values) != expected:
                raise ValueError('旧卡A/B版本与转换目标不匹配：' + name)
            base = list(ply_values)+old_values[9:]
        else:
            base = list(ply_values)+mapping['strengths']+mapping[variant+'_damage']
        p = make_props(variant, base, angle,
                       config.get('weight_plus', 0.5),
                       config.get('update', 0), config.get('fiber_cap', 1.0),
                       config.get('matrix_cap', 0.85))
        if not np.all(np.isfinite(p)) or not 0 <= p[23 if variant == 'A' else 26] <= 1:
            raise ValueError('材料包含非法参数或权重')
        new_text = '*User Material, constants='+str(len(p))
        # NLGEOM体积项及当前反馈数值切线均可能非对称。
        new_text += ', unsymm\n'
        for start in range(0, len(p), 8):
            new_text += ', '.join(format(v, '.12g') for v in p[start:start+8])+'\n'
        # 不定义自定义标签，结果只使用 SDV1、SDV2 等编号。
        depvar = f'*Depvar\n{state_count(variant)},\n'
        changes.append((constitutive.start, data_end(constitutive), depvar+new_text))
        for b in group:
            if b.key in ('depvar', 'fail stress', 'fail strain', 'user output variables'):
                changes.append((b.start, data_end(b), ''))
            if b.key in ('expansion', 'conductivity', 'specific heat'):
                raise ValueError('当前UMAT不支持热耦合材料：' + name)
        targets.add(name.lower())
        records.append({'material': name, 'ply': mapping['ply_material'], 'angle': angle,
                        'relative_elastic_error': error, 'props': p.tolist()})
    if not records:
        raise ValueError('没有找到任何目标Bin材料')
    sections = [b for b in items if b.key in ('solid section', 'shell section')
                and b.options.get('material', '').lower() in targets]
    if not sections or any(b.key != 'solid section' or 'orientation' not in b.options for b in sections):
        raise ValueError('目标材料需用于有Orientation的三维实体截面')
    # 保留原有UVARM输出请求会调用未打包的UVARM，故显式拒绝。
    for b in items:
        if b.key == 'element output' and re.search(r'\bUVARM\d*\b', b.text, re.I):
            raise ValueError('请在新导出副本中移除UVARM输出请求，改用本程序SDV')
    ordered = sorted(changes)
    if any(a[1] > b[0] for a, b in zip(ordered, ordered[1:])):
        raise AssertionError('材料替换区间重叠')
    # 切片拼接保留所有非目标字节；不要通过CAE重导出网格与取向。
    parts, position = [], 0
    for start, end, replacement in ordered:
        parts.extend([text[position:start], replacement])
        position = end
    parts.append(text[position:])
    result = select_sdv_outputs(''.join(parts), variant)
    return result, records


def convert(source, variant, config_path, output, overwrite=False):
    source, output = Path(source).resolve(), Path(output).resolve()
    if not any(output.is_relative_to(ROOT / name) for name in ('generated', 'runs')) or output == source:
        raise ValueError('输出必须位于新仓库generated/runs内，不得覆盖源文件或代码')
    reports = [output, output.with_suffix('.mapping.csv'), output.with_suffix('.report.json')]
    if not overwrite and any(p.exists() for p in reports):
        raise FileExistsError('输出已存在；请指定新名称或显式--overwrite')
    raw = source.read_bytes()
    # Abaqus关键词/数据均为ASCII；latin-1可无损保留任何既有注释字节。
    text = raw.decode('latin-1')
    config = json.loads(Path(config_path).read_text(encoding='utf-8'))
    result, records = transform(text, variant, config)
    output.parent.mkdir(parents=True, exist_ok=True)
    result = result.replace('\r\n', '\n').replace('\n', '\r\n') if b'\r\n' in raw else result
    output.write_bytes(result.encode('latin-1'))
    report = {'variant': variant, 'source': str(source),
              'source_sha256': hashlib.sha256(raw).hexdigest(),
              'output_sha256': hashlib.sha256(output.read_bytes()).hexdigest(),
              'material_count': len(records), 'state_variables': state_count(variant),
                            'output_sdv_count': len(output_indices(variant)),
              'max_relative_elastic_error': max(r['relative_elastic_error'] for r in records),
              'update': config.get('update', 0), 'materials': records}
    reports[2].write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    with reports[1].open('w', newline='', encoding='utf-8-sig') as f:
        writer = csv.writer(f)
        writer.writerow(['材料', '原始单层材料', '角度_deg', 'WCM弹性相对误差'])
        for r in records:
            writer.writerow([r['material'], r['ply'], r['angle'], r['relative_elastic_error']])
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--variant', choices=('A', 'B'), required=True)
    parser.add_argument('--config', type=Path, default=ROOT/'config/project.json')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()
    out = args.output or ROOT/'generated'/f'Job_wcm_{args.variant}.inp'
    report = convert(args.source, args.variant, args.config, out, args.overwrite)
    print(json.dumps({k: v for k, v in report.items() if k != 'materials'}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
