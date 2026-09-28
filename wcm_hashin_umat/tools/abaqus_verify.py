# -*- coding: utf-8 -*-
"""在新目录内实际运行Abaqus/Standard单元验证，保存日志和误差。

elastic：两版分别含多角度、六独立应变、倾斜Orientation及内置弹性对照。
damage：单族0/90度的纤维、基体拉压和剪切路径；检查实际SDV。
datacheck：只检查已经转换的完整气瓶输入，不进行爆破求解。
不修改系统编译配置，不在原项目写入作业文件。
"""
from pathlib import Path
import argparse
import datetime
import json
import hashlib
import shutil
import subprocess
import sys
import numpy as np
from build import build
from convert_wcm_inp import select_sdv_outputs
from model import state_count, output_indices, make_props, elastic, smeared_stiffness, rotation, tensor

ROOT = Path(__file__).resolve().parents[1]
PLY = [141000., 11400., 11400., .28, .28, .28, 7100., 7100., 7100.]
STRENGTH = [2080., 1250., 60., 290., 110., 110., 110.]
BASE = {'A': PLY+STRENGTH+[.07, .14, .2, .4, .9, .5],
        'B': PLY+STRENGTH+[133., 10., .5, 1.6, 1., .9, .5, 1e-4, 5e-3]}
CORNERS = np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],
                    [0,0,1],[1,0,1],[1,1,1],[0,1,1]], dtype=float)


def data_lines(values):
    return [' ,'.join(format(v, '.14g') for v in values[i:i+8]) for i in range(0, len(values), 8)]


def engineering(c):
    s = np.linalg.inv(c)
    return [1/s[0,0], 1/s[1,1], 1/s[2,2], -s[0,1]/s[0,0],
            -s[0,2]/s[0,0], -s[1,2]/s[1,1], c[3,3], c[4,4], c[5,5]]


def create_case(variant, suite, directory):
    lines = ['*Heading', '** WCM angle verification; mm N MPa', '*Node']
    material_lines, elements, boundary, expected = [], [], [], {}
    cases = []
    if suite == 'elastic':
        for angle in (0, 15, 37, 45, 90):
            for j in range(6):
                for builtin in (False, True):
                    cases.append((angle, j, np.eye(6)[j]*1e-5, builtin, .5))
    else:
        for angle in (0, 90):
            paths = ((0,.03),(0,-.02),(1,.012),(1,-.05)) if suite == 'completion' else ((0,.02),(0,-.012),(1,.01),(1,-.03),(3,.02),(4,.02),(5,.02))
            for j, amplitude in paths:
                direction = np.linalg.solve(elastic(PLY), np.eye(6)[j])
                local = direction/direction[j]*amplitude
                if j >= 3:
                    # 原判据在q=0处分拉/压支；有限元舍入会使纯剪切
                    # 的q正负漂移。加可忽略的拉伸偏置以固定测试分支。
                    local[1] += 1e-7
                r = rotation(angle)
                # 单族试验从给定单层应变恢复WCM方向应变。
                strain = r.T @ tensor(local, strain=True) @ r
                from model import vector
                cases.append((angle, j, vector(strain, strain=True), False, 1.))
    base_rotation = rotation(23.)
    for eid, (angle, j, e, builtin, weight) in enumerate(cases, 1):
        start = (eid-1)*8+1
        origin = np.array([2.*eid, 0., 0.])
        for k, xyz in enumerate(CORNERS*.2, start):
            point = xyz+origin
            lines.append(f'{k}, '+', '.join(map(str, point)))
            global_e = base_rotation.T @ tensor(e, strain=True) @ base_rotation
            u = global_e @ xyz
            for dof in range(3):
                boundary.append(f'{k}, {dof+1}, {dof+1}, {u[dof]:.14g}')
        elements += [f'*Element, type=C3D8, elset=E{eid}',
                     f'{eid}, '+', '.join(str(start+k) for k in range(8))]
        c = smeared_stiffness(PLY, angle, weight)
        material_lines += [f'*Material, name=M{eid}']
        if builtin:
            material_lines += ['*Elastic, type=ENGINEERING CONSTANTS']+data_lines(engineering(c))
        else:
            base = BASE[variant].copy()
            if suite == 'completion' and variant == 'B':
                base[16:20] = [4., 1.6, .045, 1.]
                base[-2:] = [.02, .02]
            p = make_props(variant, base, angle, weight=weight)
            material_lines += ['*Depvar', f'{state_count(variant)},', f'*User Material, constants={len(p)}, unsymm']+data_lines(p)
        expected[eid] = {'angle': angle, 'component': j, 'builtin': builtin,
                         'strain': list(e), 'stress_elastic': list(c @ e)}
    lines += elements
    lines += ['*Orientation, name=BASE, system=RECTANGULAR',
              ', '.join(map(str, np.r_[base_rotation[0], base_rotation[1]])), '3, 0.']
    for eid in expected:
        lines += [f'*Solid Section, elset=E{eid}, material=M{eid}, orientation=BASE', '1.,']
    lines += material_lines
    lines += ['*Step, name=LOAD, inc=1000, nlgeom=NO', '*Static',
              '0.01, 1., 1e-8, 0.01', '*Boundary']+boundary
    lines += ['*Output, field, frequency=1', '*Element Output', 'S, E, SDV',
              '*Node Output', 'U, RF', '*End Step']
    (directory/'verify.inp').write_text(select_sdv_outputs('\n'.join(lines)+'\n', variant), encoding='ascii')
    (directory/'expected.json').write_text(json.dumps(expected, indent=2), encoding='utf-8')
    return expected


def create_rotation_case(variant, directory):
    """先施加小应变，再叠加40度刚体转动，检查局部应力客观性。"""
    rbase = rotation(23.)
    f0 = rbase.T @ np.diag([1.001, 1.0004, .9999]) @ rbase
    lines = ['*Heading', '** Small strain followed by rigid rotation', '*Node']
    elems, mats, load, rotated, amplitudes = [], [], [], [], []
    for eid, angle in enumerate((15., 45.), 1):
        start = (eid-1)*8+1
        for node, x in enumerate(CORNERS*.2, start):
            lines.append(f'{node}, '+', '.join(map(str, x+np.array([eid*2.,0,0]))))
            u0 = (f0-np.eye(3)) @ x
            for j in range(3):
                load.append(f'{node}, {j+1}, {j+1}, {u0[j]:.14g}')
                name = f'A{node}_{j}'
                amplitudes.append(f'*Amplitude, name={name}')
                for t in np.linspace(0., 1., 101):
                    u = (rotation(-40.*t) @ f0-np.eye(3)) @ x
                    amplitudes.append(f'{t:.14g}, {u[j]:.14g}')
                rotated += [f'*Boundary, amplitude={name}', f'{node}, {j+1}, {j+1}, 1.']
        elems += [f'*Element, type=C3D8, elset=E{eid}',
                  f'{eid}, '+', '.join(str(start+k) for k in range(8))]
        p = make_props(variant, BASE[variant], angle)
        mats += [f'*Material, name=M{eid}', '*Depvar', f'{state_count(variant)},',
                 f'*User Material, constants={len(p)}, unsymm']+data_lines(p)
    lines += elems+['*Orientation, name=BASE',
        ', '.join(map(str, np.r_[rbase[0], rbase[1]])), '3, 0.']
    for eid in (1, 2):
        lines += [f'*Solid Section, elset=E{eid}, material=M{eid}, orientation=BASE', '1.,']
    lines += mats+amplitudes
    output = ['*Output, field, frequency=100', '*Element Output', 'S, LE, SDV', '*End Step']
    lines += ['*Step, name=LOAD, nlgeom=YES, inc=1000', '*Static',
              '0.01, 1., 1e-8, 0.01', '*Boundary']+load+output
    lines += ['*Step, name=ROTATE, nlgeom=YES, inc=1000', '*Static',
              '0.01, 1., 1e-8, 0.01']+rotated+output
    (directory/'verify.inp').write_text(select_sdv_outputs('\n'.join(lines)+'\n', variant), encoding='ascii')
    return {}


def run_command(command, directory, log, timeout=600):
    result = subprocess.run(command, cwd=directory, capture_output=True, timeout=timeout)
    (directory/log).write_bytes(result.stdout+result.stderr)
    if result.returncode or b'Traceback (most recent call last)' in result.stdout+result.stderr:
        raise RuntimeError(f'命令失败，见 {directory/log}')


def sdv(frame, element, index):
    fields = frame['fields']
    for name in (f'SDV{index}',):
        if name in fields:
            return [x['data'][0] for x in fields[name] if x['element'] == element]
    if 'SDV' in fields:
        return [x['data'][index-1] for x in fields['SDV'] if x['element'] == element]
    raise KeyError(f'ODB没有SDV{index}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--suite', choices=('elastic', 'damage', 'rotation', 'datacheck', 'completion'), default='elastic')
    parser.add_argument('--variant', choices=('A', 'B', 'both'), default='both')
    parser.add_argument('--abaqus', default=shutil.which('abaqus'))
    args = parser.parse_args()
    if not args.abaqus:
        raise RuntimeError('未找到Abaqus启动器')
    build()
    summary = []
    variants = ('A', 'B') if args.variant == 'both' else (args.variant,)
    for variant in variants:
        directory = ROOT/'runs'/('abaqus_'+args.suite+'_'+variant+'_'+datetime.datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
        directory.mkdir(parents=True)
        (directory/'scratch').mkdir()
        name = 'hashin_constant_wcm.for' if variant == 'A' else 'hashin_energy_wcm.for'
        shutil.copyfile(ROOT/'dist'/name, directory/name)
        entry = {'variant': variant, 'suite': args.suite, 'directory': str(directory),
                 'umat_sha256': hashlib.sha256((directory/name).read_bytes()).hexdigest(),
                 'passed': False}
        summary.append(entry)
        try:
            if args.suite == 'datacheck':
                shutil.copyfile(ROOT/'generated'/f'Job_wcm_{variant}.inp', directory/'verify.inp')
                expected = None
            elif args.suite == 'rotation':
                expected = create_rotation_case(variant, directory)
            else:
                expected = create_case(variant, args.suite, directory)
            command = [args.abaqus, 'job=verify', 'input=verify.inp', 'user='+name,
                       'cpus=1', 'scratch='+str(directory/'scratch'), 'interactive']
            if args.suite == 'datacheck':
                command += ['datacheck']
            run_command(command, directory, 'launcher.txt')
            logs = '\n'.join(p.read_text(errors='replace') for p in directory.glob('verify.*')
                             if p.suffix in ('.log', '.sta', '.dat'))
            launch = (directory/'launcher.txt').read_text(errors='replace')
            if 'Abaqus JOB verify COMPLETED' not in launch or 'ANALYSIS HAS NOT BEEN COMPLETED' in logs.upper():
                raise RuntimeError('Abaqus未正常完成，检查log/sta/dat')
            if args.suite != 'datacheck':
                run_command([args.abaqus, 'python', str(ROOT/'tests/extract_odb.py'),
                             'verify.odb', 'actual.json'], directory, 'extract.txt')
                frames = json.loads((directory/'actual.json').read_text())
                final = frames[-1]
                if abs(final['time']-1.) > 1e-8:
                    raise AssertionError('分析未到终点')
                if args.suite == 'elastic':
                    maximum = 0.
                    if len(final['fields']['S']) != len(expected)*8:
                        raise AssertionError('缺少元素或积分点输出')
                    for value in final['fields']['S']:
                        ref = expected[value['element']]['stress_elastic']
                        error = np.max(np.abs(np.array(value['data'])-ref))/max(1., np.max(np.abs(ref)))
                        maximum = max(maximum, float(error))
                    entry['max_elastic_error'] = maximum
                    if maximum > 1e-5:
                        raise AssertionError(f'弹性应力误差：{maximum}')
                elif args.suite == 'rotation':
                    before = {(x['element'], x['point']): np.array(x['data'])
                              for x in frames[0]['fields']['S']}
                    maximum = 0.
                    for value in final['fields']['S']:
                        ref = before[(value['element'], value['point'])]
                        error = np.max(np.abs(np.array(value['data'])-ref))/max(1., np.max(np.abs(ref)))
                        maximum = max(maximum, float(error))
                    entry['rigid_rotation_relative_error'] = maximum
                    if maximum > 1e-3:
                        raise AssertionError(f'刚体转动后局部应力变化过大：{maximum}')
                else:
                    expected_names = {f'SDV{i}' for i in output_indices(variant)}
                    actual_names = {k for k in final['fields'] if k.startswith('SDV')}
                    if actual_names != expected_names:
                        raise AssertionError(f'Unexpected ODB names: {actual_names}')
                    entry['sdv_names'] = sorted(actual_names)
                    expected_modes = [0,1,2,3]*2 if args.suite == 'completion' else [0, 1, 2, 3, 0, 0, 2]*2
                    for eid, mode in enumerate(expected_modes, 1):
                        index = 5 if mode < 2 else 6
                        if min(sdv(final, eid, index)) < (1.-1e-6 if args.suite == 'completion' else 1e-12):
                            raise AssertionError(f'Element {eid} has no damage')
                    entry['damage_checked_elements'] = len(expected)
            entry['passed'] = True
        except Exception as exc:
            entry['error'] = str(exc)
            print(entry['error'])
        (ROOT/'validation'/('abaqus_'+args.suite+'.json')).write_text(
            json.dumps(summary, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
        print(json.dumps(entry, ensure_ascii=False))
    return 0 if all(r['passed'] for r in summary) else 1


if __name__ == '__main__':
    sys.exit(main())
