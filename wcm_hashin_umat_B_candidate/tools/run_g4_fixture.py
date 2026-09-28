# -*- coding: utf-8 -*-
"""阶段4（G4）单单元 Abaqus 能量 fixture（规划9）。

- 用当前哈希 dist/hashin_energy_wcm.for 实际编译求解（不重用历史结果）；
- 单元尺寸取已证明兼容域 LC≈0.1 mm；实测 CELENT：由 SDV25（族1
  DELTA0_MT）反推 CELENT=SDV25/gamma0，不把边长直接等同 CELENT；
- NLGEOM=NO、位移控制、纯 23 剪切斜坡（eta=0、UPDATE=1，与材料点
  解析负对照同口径）；逐帧 S/LE/SDV 经 tests/extract_odb.py
  --all-frames 提取；
- 验收（规划9）：Abaqus 与同路径材料点（解析牵引恒等式）应力/状态
  归一化差 <=1e-4；能量账本 0.9*Gmt=0.45 N/mm 作为佐证（±2%）。

用法：python tools/run_g4_fixture.py --element C3D8|C3D8R|C3D6 --abaqus <bat>
"""
import argparse
import datetime
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
TOOLS = Path(__file__).resolve().parent
sys.path.insert(0, str(TOOLS))
from isolation import ensure_writable_dir, require_inside  # noqa: E402
from model import make_props  # noqa: E402

ROOT = TOOLS.parent
A = 0.1                      # 名义边长 mm（兼容域内）
GAMMA_END = 0.15             # > gamma_f = deltaf/CELENT
G23, GMT, SIGMA0 = 7100., 0.5, 110.
GAMMA0 = SIGMA0/G23
DELTA0_REF = A*GAMMA0        # 仅作初值；实测以 SDV25 为准
DELTAF = 2*GMT/SIGMA0
TARGET = 0.9*GMT
BASE25 = [141000., 11400., 11400., .28, .28, .28, 7100., 7100., 7100.,
          2080., 1250., 60., 290., 110., 110., 110.,
          133., 10., .5, 1.6, 1., .9, .5, 0., 0.]


def data_lines(values):
    return [' ,'.join(format(v, '.14g') for v in values[i:i+8])
            for i in range(0, len(values), 8)]


def build_deck(element):
    props = make_props('B', BASE25, 0, update=1)
    a = A
    if element in ('C3D8', 'C3D8R'):
        coords = [[0, 0, 0], [a, 0, 0], [a, a, 0], [0, a, 0],
                  [0, 0, a], [a, 0, a], [a, a, a], [0, a, a]]
        conn = list(range(1, 9))
        bottom, top = [1, 2, 3, 4], [5, 6, 7, 8]
    else:
        coords = [[0, 0, 0], [a, 0, 0], [0, a, 0],
                  [0, 0, a], [a, 0, a], [0, a, a]]
        conn = list(range(1, 7))
        bottom, top = [1, 2, 3], [4, 5, 6]
    lines = ['*Heading', '** WCM B G4 energy fixture; mm N MPa', '*Node']
    for k, xyz in enumerate(coords, 1):
        lines.append(f'{k}, '+', '.join(format(v, '.14g') for v in xyz))
    lines += [f'*Element, type={element}, elset=EALL',
              '1, '+', '.join(map(str, conn))]
    if element == 'C3D8R':
        lines += ['*Section Controls, name=HG_UMAT, hourglass=ENHANCED']
        lines += ['*Solid Section, elset=EALL, material=PLY_B, '
                  'controls=HG_UMAT', f'{a},']
    else:
        lines += ['*Solid Section, elset=EALL, material=PLY_B', f'{a},']
    lines += ['*Material, name=PLY_B', '*Depvar', '62,',
              '*User Material, constants=30, unsymm']+data_lines(props)
    lines += ['*Step, name=G4ENERGY, inc=10000, nlgeom=NO', '*Static',
              '0.005, 1., 1e-20, 0.1', '*Amplitude, name=RAMP',
              '0., 0.', '1., 1.', '*Boundary, amplitude=RAMP']
    for n in bottom:
        lines.append(f'{n}, 1, 6, 0.')
    for n in top:
        lines += [f'{n}, 1, 1, 0.', f'{n}, 3, 3, 0.',
                  f'{n}, 2, 2, {GAMMA_END*a:.14g}']
    lines += ['*Output, field, frequency=1', '*Element Output',
              'S, LE, SDV', '*Node Output', 'U', '*End Step']
    return '\n'.join(lines)+'\n', list(props)


def run_command(command, directory, log, timeout=900):
    result = subprocess.run(command, cwd=directory, capture_output=True,
                            timeout=timeout)
    (directory/log).write_bytes(result.stdout+result.stderr)
    if result.returncode or b'Traceback (most recent call last)' in \
            result.stdout+result.stderr:
        raise RuntimeError(f'命令失败，见 {directory/log}')


def tau_ref(gamma, lc):
    """解析牵引恒等式（材料点参考，energy套件已验证到4.5e-14）。"""
    delta = lc*gamma
    if delta < lc*GAMMA0:
        return G23*gamma
    tau_ideal = SIGMA0*(DELTAF-delta)/(DELTAF-lc*GAMMA0) if delta < DELTAF else 0.
    return 0.1*G23*gamma+0.9*max(tau_ideal, 0.)


def d_ref(gamma, lc):
    delta = lc*gamma
    d0 = lc*GAMMA0
    if delta < d0:
        return 0.
    return min(DELTAF*(delta-d0)/(delta*(DELTAF-d0)), 1.)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--element', choices=('C3D8', 'C3D8R', 'C3D6'),
                        required=True)
    parser.add_argument('--abaqus', default=r'E:\ABAQUS2025\Commands\abaqus.BAT')
    args = parser.parse_args()
    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    run_dir = ensure_writable_dir(ROOT/'runs'/f'g4_fixture_{args.element}_{stamp}')
    deck, props = build_deck(args.element)
    (run_dir/'fixture.inp').write_text(deck, encoding='ascii')
    umat = ROOT/'dist'/'hashin_energy_wcm.for'
    (run_dir/'hashin_energy_wcm.for').write_bytes(umat.read_bytes())
    command = [args.abaqus, 'job=fixture', 'input=fixture.inp',
               'user=hashin_energy_wcm.for', 'cpus=1',
               'scratch='+str(run_dir/'scratch'), 'interactive']
    (run_dir/'scratch').mkdir()
    manifest = {'element': args.element, 'gamma_end': GAMMA_END,
                'nominal_edge_mm': A, 'props': props,
                'umat_sha256': hashlib.sha256(umat.read_bytes()).hexdigest(),
                'command': command}
    (run_dir/'fixture_manifest.json').write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    with (run_dir/'launcher.txt').open('wb') as log:
        result = subprocess.run(command, cwd=run_dir, stdout=log,
                                stderr=subprocess.STDOUT, timeout=1800)
    launch = (run_dir/'launcher.txt').read_text(errors='replace')
    if result.returncode or 'Abaqus JOB fixture COMPLETED' not in launch:
        print(json.dumps({'element': args.element, 'passed': False,
                          'error': 'Abaqus作业未正常完成，见launcher.txt',
                          'run_dir': str(run_dir)}, ensure_ascii=False))
        return 1
    run_command([args.abaqus, 'python', str(ROOT/'tests'/'extract_odb.py'),
                 'fixture.odb', 'actual.json', '--all-frames', '--overwrite'],
                run_dir, 'extract.txt')
    raw = json.loads((run_dir/'actual.json').read_text())
    frames = raw['frames']
    entry = {'element': args.element, 'run_dir': str(run_dir),
             'umat_sha256': manifest['umat_sha256'], 'frames': len(frames)}
    # 逐帧取 (gamma, tau, d_sdv, delta0_sdv)：单单元单积分点。
    pts = []
    for fr in frames:
        tau = fr['fields']['S'][0]['data'][5]
        strain = fr['fields']['LE'] if 'LE' in fr['fields'] \
            else fr['fields']['E']
        gam = strain[0]['data'][5]
        d_sdv = fr['fields']['SDV8'][0]['data'][0]      # MAX_DVM（eta=0=d）
        d0_sdv = fr['fields']['SDV25'][0]['data'][0]    # 族1 DELTA0_MT
        pts.append({'time': fr['time'], 'gamma': gam, 'tau': tau,
                    'd_sdv': d_sdv, 'delta0_sdv': d0_sdv})
    onset_pts = [p for p in pts if p['delta0_sdv'] > 0]
    if not onset_pts:
        raise RuntimeError('ODB中未取到DELTA0（SDV25），无法实测CELENT')
    celement = onset_pts[0]['delta0_sdv']/GAMMA0
    entry['celem_measured_mm'] = celement
    worst_tau = worst_d = 0.
    for p in pts:
        ref = tau_ref(p['gamma'], celement)
        worst_tau = max(worst_tau, abs(p['tau']-ref)/max(1., abs(ref)))
        rd = d_ref(p['gamma'], celement)
        worst_d = max(worst_d, abs(p['d_sdv']-rd)/max(1., abs(rd)))
    W = sum(0.5*(pts[i]['tau']+pts[i+1]['tau'])*(pts[i+1]['gamma']-pts[i]['gamma'])
            for i in range(len(pts)-1))
    psi_final = 0.5*pts[-1]['tau']*pts[-1]['gamma']
    d_crack = celement*(W-psi_final)
    min_step_d = min(
        0.5*(pts[i]['tau']+pts[i+1]['tau'])*(pts[i+1]['gamma']-pts[i]['gamma'])
        -0.5*(pts[i+1]['tau']*pts[i+1]['gamma']-pts[i]['tau']*pts[i]['gamma'])
        for i in range(len(pts)-1))*celement
    entry.update({'max_tau_rel_err': worst_tau, 'max_d_rel_err': worst_d,
                  'W_density': W, 'psi_final': psi_final,
                  'D_crack_N_per_mm': d_crack,
                  'rel_vs_target_045': abs(d_crack-TARGET)/TARGET,
                  'min_step_dD': min_step_d,
                  'passed_1e-4': bool(worst_tau <= 1e-4 and worst_d <= 1e-4),
                  'passed_energy_2pct': bool(abs(d_crack-TARGET)/TARGET <= 0.02)})
    entry['passed'] = bool(entry['passed_1e-4'] and entry['passed_energy_2pct'])
    report = require_inside(ROOT/'validation'/'abaqus_candidate.json',
                            what='G4证据')
    from build import build_lock
    with build_lock(ROOT/'validation'/'.report.lock'):
        prev = json.loads(report.read_text(encoding='utf-8')) if report.is_file() else []
        prev = [e for e in prev
                if not (e.get('element') == args.element
                        and e.get('fixture', False))]
        mine = dict(entry)
        mine['fixture'] = True
        prev.append(mine)
        report.write_text(json.dumps(prev, ensure_ascii=False, indent=2)+'\n',
                          encoding='utf-8')
    print(json.dumps(entry, ensure_ascii=False, indent=2))
    return 0 if entry['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
