# -*- coding: utf-8 -*-
"""生成Abaqus诊断fixture（deck+参数清单），写入候选generated唯一目录。

规划4.3接口：--suite energy|rollback|rotation --element C3D8|C3D8R|C3D6。
本工具只生成输入与参数清单，不提交求解；实际Abaqus运行属阶段4，
由G3放行后另行执行。能量与回滚诊断不得送回prepare_test_cases的
适配流程（该流程会删除history/energy请求）。

材料卡固定为B版生产参数（规划3节），UPDATE/eta等对照改动写进
params.json并同步进deck，不静默修改。
"""
import argparse
import datetime
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from isolation import ensure_writable_dir, require_inside
from model import elastic, make_props, state_count

ROOT = Path(__file__).resolve().parents[1]
PLY = [141000., 11400., 11400., .28, .28, .28, 7100., 7100., 7100.]
STRENGTH = [2080., 1250., 60., 290., 110., 110., 110.]
BASE_B = PLY+STRENGTH+[133., 10., .5, 1.6, 1., .9, .5, 1e-4, 5e-3]
NODES_HEX = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
             [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]
NODES_WEDGE = [[0, 0, 0], [1, 0, 0], [0, 1, 0],
               [0, 0, 1], [1, 0, 1], [0, 1, 1]]


def data_lines(values):
    return [' ,'.join(format(v, '.14g') for v in values[i:i+8])
            for i in range(0, len(values), 8)]


def material_block(props, update, eta=None):
    p = np.array(props, dtype=float)
    p[-1] = update
    if eta is not None:
        p[23], p[24] = eta
    lines = ['*Depvar', f'{state_count("B")},',
             '*User Material, constants=30, unsymm']+data_lines(p)
    return lines, list(p)


def deck_header(element, length, controls):
    lines = ['*Heading', '** WCM B diagnostic fixture; mm N MPa', '*Node']
    nodes = NODES_HEX if element in ('C3D8', 'C3D8R') else NODES_WEDGE
    for k, xyz in enumerate(nodes, 1):
        point = [length*xyz[0], length*xyz[1], length*xyz[2]]
        lines.append(f'{k}, '+', '.join(map(str, point)))
    conn = list(range(1, 9)) if element in ('C3D8', 'C3D8R') else list(range(1, 7))
    lines += [f'*Element, type={element}, elset=EALL',
              '1, '+', '.join(map(str, conn))]
    lines += ['*Solid Section, elset=EALL, material=PLY_B', f'{length},']
    if controls:
        lines += ['*Section Controls, name=HG_UMAT, hourglass=ENHANCED']
        lines[-2] = ('*Solid Section, elset=EALL, material=PLY_B, '
                     'controls=HG_UMAT')
    return lines, nodes


def shear_boundary(nodes, length, gamma, step_lines):
    """纯23剪切：z=0面固定，z=L面y向位移 gamma*L，其余分量约束为零。"""
    bc = []
    for k, xyz in enumerate(nodes, 1):
        if abs(xyz[2]) < 1e-12:
            bc.append(f'{k}, 1, 6, 0.')
        else:
            bc += [f'{k}, 1, 1, 0.', f'{k}, 3, 3, 0.',
                   f'{k}, 2, 2, {gamma*length:.14g}']
    return step_lines+['*Boundary']+bc


def build_deck(suite, element):
    """返回 (deck行, 参数清单dict)。统一B版生产参数、UPDATE与eta由suite决定。"""
    length = 0.3                       # 有效长度域内（Lc_crit≈0.5868 mm）
    update = 1 if suite in ('energy', 'rotation') else 0
    eta = [0., 0.] if suite == 'energy' else [1e-4, 5e-3]
    controls = element == 'C3D8R'
    lines, nodes = deck_header(element, length, controls)
    mat, props = material_block(BASE_B, update, eta)
    lines += mat
    g23, yt = 7100., 110.
    gamma0 = yt/g23
    delta0 = length*gamma0
    sigma0 = yt
    deltaf = 2*BASE_B[18]/sigma0
    info = {'suite': suite, 'element': element, 'length_mm': length,
            'update': update, 'eta_f': eta[0], 'eta_m': eta[1],
            'props': props, 'gamma23_onset': gamma0,
            'delta0_mm': delta0, 'sigma0_mpa': sigma0,
            'deltaf_mm': deltaf,
            'gamma_end': round(1.4*deltaf/length, 8),
            'static': [0.005, 1., 1e-20, 0.1],
            'expected': ('energy: 耗能极限0.9*GMT=0.45 N/mm；'
                         'rollback: 状态不被拒绝试算污染；'
                         'rotation: 刚体转动局部应力变化<=1e-3')}
    if suite == 'energy':
        lines += ['*Step, name=SHEAR, inc=10000, nlgeom=NO', '*Static']
        lines += [', '.join(format(v, '.14g') for v in info['static'])]
        # 幅值表驱动：0→1线性升到 gamma_end，随后保持。
        lines += ['*Amplitude, name=RAMP']
        for t in np.linspace(0., 1., 41):
            lines.append(f'{t:.10g}, {t:.10g}')
        lines.append('1.1, 1.')
        bc = []
        for k, xyz in enumerate(nodes, 1):
            if abs(xyz[2]) < 1e-12:
                bc.append(f'{k}, 1, 6, 0.')
            else:
                bc += [f'{k}, 1, 1, 0.', f'{k}, 3, 3, 0.',
                       f'{k}, 2, 2, {info["gamma_end"]*length:.14g}']
        lines += ['*Boundary, amplitude=RAMP']+bc
        lines += ['*Output, field, frequency=1', '*Element Output',
                  'S, LE, SDV, SENER', '*Node Output', 'U',
                  '*Output, history, frequency=1', '*Energy Output',
                  'ALLIE, ALLSE, ALLPD, ALLAE', '*End Step']
    elif suite == 'rollback':
        lines += ['*Step, name=CYCLE, inc=10000, nlgeom=NO', '*Static']
        lines += [', '.join(format(v, '.14g') for v in info['static'])]
        lines += ['*Amplitude, name=CYC']
        for t, f in ((0., 0.), (.25, 1.), (.5, .8), (.75, 1.), (1., .9)):
            lines.append(f'{t:.10g}, {f*info["gamma_end"]:.10g}')
        bc = []
        for k, xyz in enumerate(nodes, 1):
            if abs(xyz[2]) < 1e-12:
                bc.append(f'{k}, 1, 6, 0.')
            else:
                bc += [f'{k}, 1, 1, 0.', f'{k}, 3, 3, 0.',
                       f'{k}, 2, 2, {info["gamma_end"]*length:.14g}']
        lines += ['*Boundary, amplitude=CYC']+bc
        lines += ['*Output, field, frequency=1', '*Element Output',
                  'S, LE, SDV', '*Node Output', 'U', '*End Step']
    else:  # rotation
        raise SystemExit('rotation fixture由abaqus_verify --suite rotation覆盖，'
                         '此处仅接受 energy|rollback')
    return lines, info


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--suite', choices=('energy', 'rollback'), required=True)
    parser.add_argument('--element', choices=('C3D8', 'C3D8R', 'C3D6'), required=True)
    args = parser.parse_args()
    lines, info = build_deck(args.suite, args.element)
    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    out = ensure_writable_dir(ROOT/'generated'/'diagnostics'/f'{args.suite}_{args.element}_{stamp}')
    deck = out/'fixture.inp'
    deck.write_text('\n'.join(lines)+'\n', encoding='ascii')
    params = out/'params.json'
    params.write_text(json.dumps(info, ensure_ascii=False, indent=2)+'\n',
                      encoding='utf-8')
    print(json.dumps({'deck': str(require_inside(deck)), 'params': str(params),
                      'deck_sha256': hashlib.sha256(deck.read_bytes()).hexdigest()},
                     ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
