# -*- coding: utf-8 -*-
"""阶段5-4：热点代表性材料块三档网格研究（规划10：h、h/√2、h/2）。

块体 = 单一Bin角度的均匀材料块（代表性材料块口径），沿该Bin临界IP的
全局应变状态按比例单调加载（位移场 u=t·(ε_target·x)，全部节点约束，
均匀应变精确）；物理几何/材料/载荷路径不随网格改变，仅 CELENT 随 h 变。
材料：B版生产参数（强度/G/黏性不变），UPDATE=0，C3D8R+Enhanced。

对比（规划10）：起始时刻、响应曲线、归一化耗能——细两档差<=5%；
C3D8R 的 ALLAE/ALLIE<=5%。三档后仍不收敛/报错 => BLOCKED_MESH 证据。

用法：python tools/run_mesh_block_study.py --bin 40 --abaqus <bat>
"""
import argparse
import csv
import datetime
import hashlib
import json
import math
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
RUN = ROOT/'runs'/'mesh_survey_20260927'
BLOCK_MM = 2.0
BASE25 = [141000., 11400., 11400., .28, .28, .28, 7100., 7100., 7100.,
          2080., 1250., 60., 290., 110., 110., 110.,
          133., 10., .5, 1.6, 1., .9, .5, 1e-4, 5e-3]


def data_lines(values):
    return [' ,'.join(format(v, '.14g') for v in values[i:i+8])
            for i in range(0, len(values), 8)]


def pick_median_point(bin_id, mode='MC'):
    rows = [r for r in csv.DictReader(open(RUN/'rho_critical.csv',
                                           encoding='utf-8'))
            if r['mode'] == mode and r['bin'] == str(bin_id)]
    if not rows:
        raise RuntimeError(f'Bin {bin_id} 无 {mode} 临界点')
    rows.sort(key=lambda r: float(r['rho']))
    med = rows[len(rows)//2]
    return int(med['element']), float(med['rho']), med['family']


def element_strain(element, time_key='last_converged'):
    raw = json.loads((RUN/'strain_field.json').read_text())
    frames = [f for f in raw['frames'] if f['ips']]
    frames.sort(key=lambda f: f['time'])
    if time_key == 'last_converged':
        fr = frames[-2] if len(frames) >= 2 and abs(
            frames[-1]['time']-frames[-2]['time']) < 1e-9 else frames[-1]
    else:
        fr = frames[-1]
    for row in fr['ips']:
        if row[0] == element:
            return np.array(row[1:7]), fr['time']
    raise RuntimeError(f'单元 {element} 不在应变场中')


def build_deck(h, theta, eps_target):
    props = make_props('B', BASE25, theta, update=0)
    n = max(2, int(round(BLOCK_MM/h)))
    h_eff = BLOCK_MM/n
    lines = ['*Heading', '** WCM B mesh block study; mm N MPa', '*Node']
    nid = {}
    k = 0
    for iz in range(n+1):
        for iy in range(n+1):
            for ix in range(n+1):
                k += 1
                nid[(ix, iy, iz)] = k
                lines.append(f'{k}, {ix*h_eff:.14g}, {iy*h_eff:.14g}, '
                             f'{iz*h_eff:.14g}')
    lines += [f'*Element, type=C3D8R, elset=EALL']
    conn_lines = []
    for iz in range(n):
        for iy in range(n):
            for ix in range(n):
                c = (nid[(ix, iy, iz)], nid[(ix+1, iy, iz)],
                     nid[(ix+1, iy+1, iz)], nid[(ix, iy+1, iz)],
                     nid[(ix, iy, iz+1)], nid[(ix+1, iy, iz+1)],
                     nid[(ix+1, iy+1, iz+1)], nid[(ix, iy+1, iz+1)])
                conn_lines.append(f'{len(conn_lines)+1}, '
                                  + ', '.join(map(str, c)))
    lines += conn_lines
    lines += ['*Section Controls, name=HG_UMAT, hourglass=ENHANCED',
              '*Solid Section, elset=EALL, material=PLY_B, '
              'controls=HG_UMAT', f'{h_eff},',
              '*Material, name=PLY_B', '*Depvar', '62,',
              '*User Material, constants=30, unsymm']+data_lines(props)
    lines += ['*Step, name=MESH, inc=10000, nlgeom=NO', '*Static',
              '0.005, 1., 1e-20, 0.1', '*Amplitude, name=RAMP',
              '0., 0.', '1., 1.', '*Boundary, amplitude=RAMP']
    # 均匀应变：u_i = Σ_j ε_ij·x_j（工程剪应变，对称张量），全部节点约束
    e11, e22, e33, e12, e13, e23 = eps_target
    for (ix, iy, iz), label in nid.items():
        x, y, z = ix*h_eff, iy*h_eff, iz*h_eff
        u = (e11*x+e12*y+e13*z, e12*x+e22*y+e23*z, e13*x+e23*y+e33*z)
        for dof in range(3):
            lines.append(f'{label}, {dof+1}, {dof+1}, {u[dof]:.14g}')
    lines += ['*Output, field, frequency=1', '*Element Output',
              'S, E, SDV', '*Node Output', 'U',
              '*Output, history, frequency=1', '*Energy Output',
              'ALLSE, ALLIE, ALLPD, ALLAE', '*End Step']
    return '\n'.join(lines)+'\n', h_eff, n


def run_command(command, directory, log, timeout=1800):
    result = subprocess.run(command, cwd=directory, capture_output=True,
                            timeout=timeout)
    (directory/log).write_bytes(result.stdout+result.stderr)
    if result.returncode or b'Traceback (most recent call last)' in \
            result.stdout+result.stderr:
        raise RuntimeError(f'命令失败，见 {directory/log}')


def extract(directory, abaqus):
    run_command([abaqus, 'python', str(ROOT/'tests'/'extract_odb.py'),
                 'block.odb', 'actual.json', '--all-frames', '--overwrite'],
                directory, 'extract.txt')
    raw = json.loads((directory/'actual.json').read_text())
    frames = raw['frames']
    run_command([abaqus, 'python', str(ROOT/'tests'
                 /'extract_history_energy.py'),
                 'block.odb', 'history.json'], directory, 'extract_h.txt')
    hist = json.loads((directory/'history.json').read_text())
    return frames, hist


def ledger(frames):
    """逐帧账本（代表性单元1）：W=∫σ:dε，Ψ=0.5σ·ε，D=W-Ψ（体积密度）。"""
    pts = []
    for fr in frames:
        tau = np.array(fr['fields']['S'][0]['data'])
        eps = np.array((fr['fields'].get('E')
                        or fr['fields'].get('LE'))[0]['data'])
        d_mt = fr['fields']['SDV13'][0]['data'][0]
        d_mc = fr['fields']['SDV14'][0]['data'][0]
        init = [fr['fields'][f'SDV{i}'][0]['data'][0] for i in (19, 20, 21, 22)]
        pts.append({'t': fr['time'], 'tau': tau, 'eps': eps,
                    'd_mt': d_mt, 'd_mc': d_mc, 'init': init})
    W = 0.
    for i in range(len(pts)-1):
        W += 0.5*float(np.dot(pts[i]['tau']+pts[i+1]['tau'],
                              pts[i+1]['eps']-pts[i]['eps']))
    psi_final = 0.5*float(np.dot(pts[-1]['tau'], pts[-1]['eps']))
    return {'W_density': W, 'psi_final': psi_final,
            'D_density': W-psi_final, 'points': pts}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bin', type=int, default=40)
    parser.add_argument('--mode', default='MC')
    parser.add_argument('--levels', default='0.5,0.3536,0.25')
    parser.add_argument('--block', type=float, default=2.0)
    parser.add_argument('--abaqus', default=r'E:\ABAQUS2025\Commands\abaqus.BAT')
    args = parser.parse_args()
    global BLOCK_MM
    BLOCK_MM = args.block
    element, rho_survey, family = pick_median_point(args.bin, args.mode)
    eps_target, t_state = element_strain(element)
    rows = [r for r in csv.DictReader(open(RUN/'rho_critical.csv',
                                           encoding='utf-8'))
            if int(r['element']) == element]
    theta = float(rows[0]['angle_deg'])
    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    entries = []
    for h in [float(x) for x in args.levels.split(',')]:
        run_dir = ensure_writable_dir(ROOT/'runs'/f'mesh_block_{args.bin}_'
                                      f'{h:.4f}_{stamp}'.replace('.', 'p'))
        deck, h_eff, n = build_deck(h, theta, eps_target)
        (run_dir/'block.inp').write_text(deck, encoding='ascii')
        umat = ROOT/'dist'/'hashin_energy_wcm.for'
        (run_dir/'hashin_energy_wcm.for').write_bytes(umat.read_bytes())
        command = [args.abaqus, 'job=block', 'input=block.inp',
                   'user=hashin_energy_wcm.for', 'cpus=1',
                   'scratch='+str(run_dir/'scratch'), 'interactive']
        (run_dir/'scratch').mkdir()
        (run_dir/'block_manifest.json').write_text(json.dumps(
            {'h_requested': h, 'h_effective': h_eff, 'n_elements_per_side': n,
             'bin': args.bin, 'mode': args.mode, 'theta': theta,
             'element_survey': element, 'rho_survey': rho_survey,
             'eps_target': list(eps_target),
             'umat_sha256': hashlib.sha256(umat.read_bytes()).hexdigest(),
             'command': command}, ensure_ascii=False, indent=2)+'\n',
            encoding='utf-8')
        with (run_dir/'launcher.txt').open('wb') as log:
            result = subprocess.run(command, cwd=run_dir, stdout=log,
                                    stderr=subprocess.STDOUT, timeout=3600)
        launch = (run_dir/'launcher.txt').read_text(errors='replace')
        completed = 'Abaqus JOB block COMPLETED' in launch \
            and result.returncode == 0
        entry = {'h': h, 'h_effective': h_eff, 'n': n,
                 'run_dir': str(run_dir), 'completed': completed}
        if completed:
            frames, hist = extract(run_dir, args.abaqus)
            led = ledger(frames)
            entry['frames'] = len(frames)
            entry['W_density'] = led['W_density']
            entry['D_density'] = led['D_density']
            entry['t_onset'] = next((p['t'] for p in led['points']
                                     if any(v >= 1 for v in p['init'])), None)
            entry['d_mt_final'] = led['points'][-1]['d_mt']
            entry['d_mc_final'] = led['points'][-1]['d_mc']
            ae = hist.get('ALLAE', [])
            ie = hist.get('ALLIE', [])
            if ae and ie:
                ratios = [abs(a)/(abs(i)+1e-30) for a, i in zip(ae, ie)
                          if abs(i) > 1e-9]
                entry['ALLAE_over_ALLIE_max'] = max(ratios) if ratios else 0.0
        entries.append(entry)
        print(json.dumps(entry, ensure_ascii=False))
    # 细两档对比（规划10：差<=5%）
    ok = len(entries) >= 2 and entries[-1].get('completed') \
        and entries[-2].get('completed')
    comparison = {'levels_compared': None, 'passed_5pct': None}
    if ok:
        a, b = entries[-2], entries[-1]
        diffs = {
            'onset_time': abs((a['t_onset'] or 0)-(b['t_onset'] or 0))
            / max(a['t_onset'] or 1e-30, 1e-30),
            'W_density': abs(a['W_density']-b['W_density'])
            / max(abs(b['W_density']), 1e-30),
            'D_density': abs(a['D_density']-b['D_density'])
            / max(abs(b['D_density']), 1e-30)}
        comparison = {'levels_compared': [a['h'], b['h']],
                      'diffs': diffs,
                      'passed_5pct': all(v <= 0.05 for v in diffs.values()),
                      'passed_allae': (entries[-1].get('ALLAE_over_ALLIE_max')
                                       or 0) <= 0.05}
    report = {'bin': args.bin, 'mode': args.mode, 'theta': theta,
              'element_survey': element, 'rho_survey': rho_survey,
              'eps_target': list(eps_target), 'levels': entries,
              'comparison': comparison,
              'verdict': ('PASS' if comparison.get('passed_5pct')
                          and comparison.get('passed_allae') else
                          'BLOCKED_MESH' if not ok else 'FAIL')}
    out = require_inside(ROOT/'validation'/'mesh_block_study.json')
    out.write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n',
                   encoding='utf-8')
    print(json.dumps({'verdict': report['verdict'],
                      'comparison': comparison}, ensure_ascii=False, indent=1))
    return 0 if report['verdict'] == 'PASS' else 1


if __name__ == '__main__':
    sys.exit(main())
