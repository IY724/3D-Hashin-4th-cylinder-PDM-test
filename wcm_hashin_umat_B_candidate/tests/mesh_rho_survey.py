# -*- coding: utf-8 -*-
"""阶段5-3：ρ=δf/δ0 全场普查（规划10：逐族/逐模式，不用单一长度要求全瓶）。

输入：
  runs/mesh_survey_*/strain_field.json  （fallback ODB 提取的 LE 应变场，
        NLGEOM=NO 下与工程应变差 O(ε²)≈5e-5，对普查可忽略）
  geometry_zhuning_joint_rebuild/WCM_Job.inp  （节点/单元/Bin映射，只读）
  geometry_zhuning_joint_rebuild/WCM_Zhuning_Joint_Rebuilt_WindAngles.ang
输出：
  validation/mesh_rho_survey.json       汇总
  runs/mesh_survey_*/rho_critical.csv   ρ<=1.2 的关键点清单

口径（与本构账本一致，Python独立实现）：
  ε_L = B(±θ)·ε_G；EFF=C0·ε_L；DEL/SIG 按 WCM_EQ；FI 按 WCM_HASHIN；
  ρ_I = δf_I/δ0_I = (2·G_I/SIG_I)/(LC·|等效位移|)，仅 DEL>0 且 SIG>0 的点。
  LC 取每单元体积^(1/3)（G4 fixture 实测 CELENT=体积^(1/3)），与 fallback
  日志的 L 逐元交叉核对。LC_req = 2·G_I/(1.2·SIG_I·γ_eq)（ρ>=1.2 工程裕量）。
"""
import csv
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

TOOLS = Path(__file__).resolve().parent.parent/'tools'
sys.path.insert(0, str(TOOLS))
from model import elastic  # noqa: E402
from prepare_test_cases import tank_topology, read_angles  # noqa: E402

ROOT = TOOLS.parent
REPO = ROOT.parent
RUN = ROOT/'runs'/'mesh_survey_20260927'
G = {'FT': 133.0, 'FC': 10.0, 'GMT': 0.5, 'GMC': 1.6}
C0 = elastic([141000., 11400., 11400., .28, .28, .28, 7100., 7100., 7100.])
STRENGTH = {'XT': 2080., 'XC': 1250., 'YT': 60., 'YC': 290.,
            'S12': 110., 'S23': 110.}
RHO_TARGET = 1.2


def parse_inp_nodes_elements():
    """从 rebuild INP 解析 TANK-1 的节点坐标与单元连接（只读）。"""
    text = (REPO/'geometry_zhuning_joint_rebuild'/'WCM_Job.inp'
            ).read_bytes().decode('latin-1')
    lines = text.splitlines()
    nodes, elems, current = {}, {}, None
    for line in lines:
        head = line.strip().lower()
        if head.startswith('*node'):
            current = 'n'
            continue
        if head.startswith('*element, type=c3d'):
            current = 'e'
            continue
        if head.startswith('*'):
            current = None
            continue
        if current == 'n' and ',' in line:
            parts = [x.strip() for x in line.split(',') if x.strip()]
            if len(parts) == 4:
                nodes[int(parts[0])] = np.array([float(parts[1]),
                                                 float(parts[2]),
                                                 float(parts[3])])
        elif current == 'e' and ',' in line:
            parts = [int(x.strip()) for x in line.split(',') if x.strip()]
            if len(parts) >= 7:
                elems[parts[0]] = parts[1:]
    return nodes, elems


def element_volume(nodes, conn):
    """C3D8: 绕对角0-6的6四面体；C3D6: 3四面体。"""
    pts = [nodes[c] for c in conn]
    if len(conn) == 8:
        order = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6),
                 (0, 6, 7, 5), (0, 6, 5, 4), (0, 6, 4, 1)]
    else:
        order = [(0, 1, 2, 4), (0, 2, 5, 4), (0, 4, 5, 3)]
    vol = 0.0
    for a, b, c, d in order:
        mat = np.array([pts[b]-pts[a], pts[c]-pts[a], pts[d]-pts[a]])
        vol += abs(np.linalg.det(mat))/6.0
    return vol


def rot6(theta):
    t = math.radians(theta)
    c, s = math.cos(t), math.sin(t)
    if abs(c) < 1e-14:
        c = 0.0
    if abs(s) < 1e-14:
        s = 0.0
    return np.array([
        [c*c, s*s, 0, c*s, 0, 0],
        [s*s, c*c, 0, -c*s, 0, 0],
        [0, 0, 1, 0, 0, 0],
        [-2*c*s, 2*c*s, 0, c*c-s*s, 0, 0],
        [0, 0, 0, 0, c, s],
        [0, 0, 0, 0, -s, c]])


def batch_rot6(thetas):
    t = np.radians(thetas)
    c, s = np.cos(t), np.sin(t)
    c = np.where(np.abs(c) < 1e-14, 0.0, c)
    s = np.where(np.abs(s) < 1e-14, 0.0, s)
    n = len(thetas)
    b = np.zeros((n, 6, 6))
    b[:, 0, 0] = c*c
    b[:, 0, 1] = s*s
    b[:, 0, 3] = c*s
    b[:, 1, 0] = s*s
    b[:, 1, 1] = c*c
    b[:, 1, 3] = -c*s
    b[:, 2, 2] = 1.0
    b[:, 3, 0] = -2*c*s
    b[:, 3, 1] = 2*c*s
    b[:, 3, 3] = c*c-s*s
    b[:, 4, 4] = c
    b[:, 4, 5] = s
    b[:, 5, 4] = -s
    b[:, 5, 5] = c
    return b


def eq_del_sig(e, s):
    """WCM_EQ 的向量化独立实现。e,s: (N,6)。返回 DEL,SIG: (N,4)。"""
    sh = e[:, 3]**2+e[:, 4]**2+e[:, 5]**2
    shw = s[:, 3]*e[:, 3]+s[:, 4]*e[:, 4]+s[:, 5]*e[:, 5]
    d = np.stack([
        np.sqrt(np.maximum(e[:, 0], 0)**2+e[:, 3]**2+e[:, 4]**2),
        np.maximum(-e[:, 0], 0),
        np.sqrt(np.maximum(e[:, 1], 0)**2+np.maximum(e[:, 2], 0)**2+sh),
        np.sqrt(np.maximum(-e[:, 1], 0)**2+np.maximum(-e[:, 2], 0)**2+sh)],
        axis=1)
    with np.errstate(divide='ignore', invalid='ignore'):
        g0 = np.where(d[:, 0] > 0,
                      (np.maximum(s[:, 0], 0)*np.maximum(e[:, 0], 0)
                       + s[:, 3]*e[:, 3]+s[:, 4]*e[:, 4])
                      / np.where(d[:, 0] > 0, d[:, 0], 1), 0)
        g1 = np.where(d[:, 1] > 0,
                      np.maximum(-s[:, 0], 0)*np.maximum(-e[:, 0], 0)
                      / np.where(d[:, 1] > 0, d[:, 1], 1), 0)
        g2 = np.where(d[:, 2] > 0,
                      (np.maximum(s[:, 1], 0)*np.maximum(e[:, 1], 0)
                       + np.maximum(s[:, 2], 0)*np.maximum(e[:, 2], 0)+shw)
                      / np.where(d[:, 2] > 0, d[:, 2], 1), 0)
        g3 = np.where(d[:, 3] > 0,
                      (np.maximum(-s[:, 1], 0)*np.maximum(-e[:, 1], 0)
                       + np.maximum(-s[:, 2], 0)*np.maximum(-e[:, 2], 0)+shw)
                      / np.where(d[:, 3] > 0, d[:, 3], 1), 0)
    return d, np.stack([g0, g1, g2, g3], axis=1)


def fi_branch(s):
    """分支感知 FI（WCM_HASHIN 同口径）。s: (N,6)。返回 (N,4)。"""
    q = s[:, 1]+s[:, 2]
    sh = (s[:, 5]**2-s[:, 1]*s[:, 2])/STRENGTH['S23']**2 \
        + (s[:, 3]/STRENGTH['S12'])**2 + (s[:, 4]/STRENGTH['S12'])**2
    fi = np.zeros((len(s), 4))
    pos = s[:, 0] >= 0
    fi[pos, 0] = (s[pos, 0]/STRENGTH['XT'])**2+sh[pos]
    fi[~pos, 1] = (s[~pos, 0]/STRENGTH['XC'])**2
    qp = q >= 0
    fi[qp, 2] = (q[qp]/STRENGTH['YT'])**2+sh[qp]
    fi[~qp, 3] = ((STRENGTH['YC']/(2*STRENGTH['S23']))**2-1) \
        * q[~qp]/STRENGTH['YC'] + (q[~qp]/(2*STRENGTH['S23']))**2+sh[~qp]
    return fi


def parse_fallback_lc():
    """从 fallback msg 解析 element->L（Abaqus自己的CELENT），供交叉核对。"""
    msg = (REPO/'test'/'B_energy_softening'/'runs'
           /'B_20260927_153131_fallback'/'Rebuild_PDM_B.msg')
    if not msg.is_file():
        return {}
    text = msg.read_text(encoding='utf-8', errors='replace')
    pat = re.compile(
        r'WCM_TANK1_MAT1_BIN(\d+)\s*\n\s*([-\d.E+]+)\s+[-\d.E+]+\s+[-\d.E+]+\s*\n'
        r'\s*WCM_FALLBACK element/point/family/mode:\s*(\d+)')
    lc = {}
    for m in pat.finditer(text):
        elem = int(m.group(3))
        lc.setdefault(elem, float(m.group(2)))
    return lc


def main():
    run = RUN
    field = json.loads((run/'strain_field.json').read_text())
    # --- 几何与映射 ---
    nodes, elems = parse_inp_nodes_elements()
    deck = (REPO/'geometry_zhuning_joint_rebuild'/'WCM_Job.inp'
            ).read_bytes().decode('latin-1')
    _, sections, _, _, resolve = tank_topology(deck)
    comp_elements = {}
    for sec in sections:
        mat = sec.options['material']
        m = re.match(r'WCM_Tank1_Mat1_Bin(\d+)', mat, re.I)
        for label in resolve(sec.options['elset']):
            comp_elements[label] = int(m.group(1))
    ang = read_angles(REPO/'geometry_zhuning_joint_rebuild'
                      /'WCM_Zhuning_Joint_Rebuilt_WindAngles.ang',
                      set(elems))
    log_lc = parse_fallback_lc()
    # --- 每单元 LC：优先 fallback 日志的 Abaqus CELENT，否则体积^(1/3)；交叉核对 ---
    lc_vol, lc_err, lc_src = {}, [], {}
    for label, conn in elems.items():
        if label in comp_elements:
            v = element_volume(nodes, conn)
            lc_vol[label] = v**(1/3)
            lc_src[label] = 'vol'
            if label in log_lc and log_lc[label] > 0:
                rel = abs(lc_vol[label]-log_lc[label])/log_lc[label]
                lc_err.append(rel)
                if rel > 0.2:
                    lc_src[label] = 'log'
                    lc_vol[label] = log_lc[label]
    lc_check = {'pairs': len(lc_err),
                'max_rel_diff': max(lc_err) if lc_err else None,
                'pairs_over_20pct': sum(1 for x in lc_err if x > 0.2),
                'rule': '体积^(1/3)与fallback日志CELENT差>20%的单元改用日志值',
                'note': 'CELENT=体积^(1/3)已经G4 fixture实测（立方/楔形体）'}
    # --- 逐帧普查（同时间双帧：首帧=收敛帧，次帧=中止尝试态，标注区分）---
    modes = ('FT', 'FC', 'MT', 'MC')
    summary = {'frames': [], 'lc_check': lc_check,
               'composite_elements': len(comp_elements),
               'note': 'LE应变（NLGEOM=NO下与工程应变差O(e^2)）；'
                       '状态为fallback分支(即时完成)收敛解，区域/尺度为近似。'}
    critical_rows = []
    seen_time = {}
    for fr in field['frames']:
        t = fr['time']
        seen_time[t] = seen_time.get(t, 0)+1
        frame_tag = 'converged' if seen_time[t] == 1 else 'unconverged_attempt'
        ips = {}
        for row in fr['ips']:
            ips.setdefault(row[0], []).append(row[1:7])
        # 组装向量
        labels, eps_g, bins, thetas, lcs = [], [], [], [], []
        for elem, strains in ips.items():
            if elem not in comp_elements or elem not in lc_vol:
                continue
            for e6 in strains:
                labels.append(elem)
                eps_g.append(e6)
                bins.append(comp_elements[elem])
                thetas.append(math.degrees(ang[elem]))
                lcs.append(lc_vol[elem])
        eps_g = np.array(eps_g)
        lcs = np.array(lcs)
        bins = np.array(bins)
        labels = np.array(labels)
        frame_stat = {'time': t, 'tag': frame_tag, 'points': len(eps_g),
                      'modes': {}}
        crit = []
        for fam, sign in (('PLUS', 1.0), ('MINUS', -1.0)):
            thetas_f = np.array(thetas)*sign
            bl = batch_rot6(thetas_f)
            eps_l = np.einsum('nij,nj->ni', bl, eps_g)
            eff = np.einsum('ij,nj->ni', C0, eps_l)
            d, g = eq_del_sig(eps_l, eff)
            fi = fi_branch(eff)
            for mi, mname in enumerate(modes):
                gv = [G['FT'], G['FC'], G['GMT'], G['GMC']][mi]
                valid = (d[:, mi] > 0) & (g[:, mi] > 0)
                rho = np.full(len(d), np.nan)
                rho[valid] = (2*gv/g[valid, mi])/d[valid, mi]
                gamma_eq = d[valid, mi]/lcs[valid]
                lc_req = np.full(len(d), np.nan)
                lc_req[valid] = 2*gv/(RHO_TARGET*g[valid, mi]*gamma_eq)
                over = fi[:, mi] >= 1
                n_valid = int(valid.sum())
                n_over = int((valid & over).sum())
                n_le1 = int((valid & (rho <= 1)).sum())
                n_le12 = int((valid & (rho <= RHO_TARGET)).sum())
                rho_over = rho[valid & over]
                stat = {'family': fam, 'mode': mname,
                        'points_del_gt0': n_valid,
                        'points_fi_ge1': n_over,
                        'points_rho_le1': n_le1,
                        'points_rho_le1.2': n_le12}
                if len(rho_over):
                    stat['rho_over_percentiles'] = {
                        'p0': float(np.min(rho_over)),
                        'p10': float(np.percentile(rho_over, 10)),
                        'p50': float(np.percentile(rho_over, 50)),
                        'p90': float(np.percentile(rho_over, 90))}
                # ρ<=1（不相容）子集的全总体统计
                incomp = valid & (rho <= 1)
                if incomp.sum():
                    req1 = rho[incomp]*lcs[incomp]/1.0
                    stat['lc_req_rho_le1_mm'] = {
                        'n': int(incomp.sum()),
                        'p50': float(np.percentile(req1, 50)),
                        'p90': float(np.percentile(req1, 90)),
                        'min': float(np.min(req1)),
                        'max': float(np.max(req1))}
                # ρ<=1.2 点的细化尺度（可行动口径）：LC_req = LC*rho/1.2
                band = valid & (rho <= RHO_TARGET)
                if band.sum():
                    req = rho[band]*lcs[band]/RHO_TARGET
                    stat['lc_req_rho_le1.2_mm'] = {
                        'p50': float(np.percentile(req, 50)),
                        'p90': float(np.percentile(req, 90)),
                        'max': float(np.max(req)),
                        'min': float(np.min(req))}
                    stat['lc_req_bins_top'] = {
                        str(int(b)): float(np.percentile(
                            rho[band & (bins == b)]*lcs[band & (bins == b)]
                            /RHO_TARGET, 90))
                        for b in np.unique(bins[band])}
                frame_stat['modes'][fam+'_'+mname] = stat
                if n_le12:
                    idx = np.where(valid & (rho <= RHO_TARGET))[0]
                    order = idx[np.argsort(rho[idx])]
                    for i in order[:2000]:
                        crit.append({
                            'time': t, 'tag': frame_tag, 'family': fam,
                            'mode': mname,
                            'element': int(labels[i]), 'bin': int(bins[i]),
                            'angle_deg': round(thetas[i], 4),
                            'lc_mm': round(float(lcs[i]), 4),
                            'rho': round(float(rho[i]), 4),
                            'lc_req_mm': round(float(lc_req[i]), 4),
                            'fi_ge1': bool(over[i])})
        summary['frames'].append(frame_stat)
        if frame_tag == 'converged':
            critical_rows.extend(crit)
    # --- 已知9单元交叉核对 ---
    known = [13426, 13499, 13572, 50964, 51037, 51110, 13574, 13428, 13501]
    known_rows = [r for r in critical_rows if r['element'] in known]
    summary['known_elements'] = {
        'elements': known,
        'rows_found': len(known_rows),
        'note': 'fallback日志取证值 d0=0.0152-0.0296/d f=0.0094-0.0140 对应 '
                'rho=0.44-0.76（LC≈2.3mm，增量端点口径）；本普查为当前状态口径。'}
    (run/'rho_critical.csv').write_text(
        'time,family,mode,element,bin,angle_deg,lc_mm,rho,lc_req_mm,fi_ge1\n'
        + '\n'.join(','.join(str(r[k]) for k in
                             ('time', 'family', 'mode', 'element', 'bin',
                              'angle_deg', 'lc_mm', 'rho', 'lc_req_mm',
                              'fi_ge1')) for r in critical_rows)+'\n',
        encoding='utf-8')
    out = ROOT/'validation'/'mesh_rho_survey.json'
    out.write_text(json.dumps(summary, ensure_ascii=False, indent=1)+'\n',
                   encoding='utf-8')
    print(json.dumps({'frames': len(summary['frames']),
                      'critical_rows': len(critical_rows),
                      'lc_check': lc_check}, ensure_ascii=False, indent=1))


if __name__ == '__main__':
    main()
