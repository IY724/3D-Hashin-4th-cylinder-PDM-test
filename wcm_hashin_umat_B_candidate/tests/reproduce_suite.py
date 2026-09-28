# -*- coding: utf-8 -*-
"""阶段1复现套件：重编译严格版负对照 + 旧版静默跳过对照 + 解析基准。

规划5.1：
- 用候选源码重新编译（不用历史DLL），在 0°、纯剪 gamma12=0.006/0.01/0.0155、
  Lc=2.23 mm 路径上复现严格版 1033（族1基体拉伸 deltaf<=delta0）。
- 旧程序（baseline/hashin.for）同路径静默跳过：损伤冻结、不报错，仅作对照。
- 解析基准：gamma0=110/7100；该路径 Lc_crit=2*7100*0.5/110^2。
- 端点起始的步长依赖证据：同终点的 1 步与多步加载得到不同起始量，
  即阶段2要修的缺陷。
输出：validation/reproduction.json + 运行目录逐接受点 CSV。
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

TOOLS = Path(__file__).resolve().parent
sys.path.insert(0, str(TOOLS))
from model import make_props  # noqa: E402

PLY = [141000., 11400., 11400., .28, .28, .28, 7100., 7100., 7100.]
STRENGTH = [2080., 1250., 60., 290., 110., 110., 110.]
BASE_B = PLY+STRENGTH+[133., 10., .5, 1.6, 1., .9, .5, 1e-4, 5e-3]
XT, XC, YT, YC, S12, S13, S23 = STRENGTH   # PROPS(10:16)槽序
G12, G23 = 7100., 7100.
GMT = 0.5
GAMMA0_23 = S23/G23               # 23纯剪基体拉伸起始：0.0154929577...
GAMMA0_12 = S12/G12               # 12纯剪起始（S12=110，与23同为110/7100）
LC_CRIT_MT_SHEAR = 2*G23*GMT/S23**2  # 0.5867768595...（该路径条件，非全瓶建议）


def ramp(core, p, gamma_target, lc, increments, shear=3):
    """0°单层纯剪增量加载；逐接受点记录。shear=3 表示 gamma12(E4)。"""
    state = np.zeros(84)
    e0 = np.zeros(6)
    rows = []
    error = None
    de = np.zeros(6)
    de[shear] = gamma_target/increments
    for i in range(1, increments+1):
        try:
            state, s, c, err = core.update('B', p, state, e0, de,
                                           length=lc, allow_error=True)
        except RuntimeError:
            raise
        row = {'increment': i, 'gamma': float(de[shear]*i), 'ierr': err,
               'FI_FT': float(state[24]), 'FI_FC': float(state[25]),
               'FI_MT': float(state[26]), 'FI_MC': float(state[27]),
               'delta0_FT': float(state[12]), 'delta0_MT': float(state[14]),
               'sigma0_FT': float(state[16]), 'sigma0_MT': float(state[18]),
               'd_FT': float(state[0]), 'd_MT': float(state[2]),
               'dv_MT': float(state[22]),
               'tau': float(s[shear])}
        rows.append(row)
        if err:
            error = err
            break
        e0 = e0+de
    return rows, error


def single_increment(core, p, gamma, lc, shear=3):
    state = np.zeros(84)
    de = np.zeros(6)
    de[shear] = gamma
    state, s, c, err = core.update('B', p, state, np.zeros(6), de,
                                   length=lc, allow_error=True)
    return state, s, err


def legacy_control(core, gamma, lc, shear=3):
    """旧程序同路径：应静默跳过（不报错、MT损伤冻结为0）。"""
    old = np.zeros(24)
    e0 = np.zeros(6)
    de = np.zeros(6)
    de[shear] = gamma
    new, s, c = core.original('B', old, e0, de, length=lc)
    return {'legacy_d_mt': float(new[2]), 'legacy_init_mt': float(new[10]),
            'legacy_delta0_mt': float(new[14]), 'legacy_sigma0_mt': float(new[18]),
            'legacy_kappa_mt': float(new[6]), 'legacy_tau': float(s[shear]),
            'silent_skip': bool(new[10] >= 1 and new[2] == 0.0)}


def strict_ply_forensics(core, gamma, lc, shear_index):
    """直接调用候选版 WCM_PLY 取族内局部状态（报错路径上 WCM_EVAL 不回拷
    LOCAL，族内 delta0/sigma0 只有从 WCM_PLY 才可见）。
    候选签名：WCM_PLY(MODE,P,C0,OLD,E0,E,DT,LC,IADV,IONSET,NEW,S,C,IERR)。"""
    import ctypes as ct
    from model import elastic as ply_elastic
    c0 = np.asfortranarray(ply_elastic(PLY))
    e = np.zeros(6)
    e[shear_index] = gamma
    e0 = np.zeros(6)
    old = np.zeros(40)
    new = np.zeros(40)
    s = np.zeros(6)
    c = np.zeros((6, 6), order='F')
    p = np.array(make_props('B', BASE_B, 0), dtype=np.float64)
    mode, iadv, ionset, ierr = (ct.c_int(2), ct.c_int(1), ct.c_int(1),
                                ct.c_int(0))
    dt, lc_ = ct.c_double(0.01), ct.c_double(lc)

    def ptr(a):
        return a.ctypes.data_as(ct.POINTER(ct.c_double))

    core.library.wcm_ply_(ct.byref(mode), ptr(p), ptr(c0),
                          ptr(old), ptr(e0), ptr(e), ct.byref(dt), ct.byref(lc_),
                          ct.byref(iadv), ct.byref(ionset), ptr(new), ptr(s),
                          ptr(c), ct.byref(ierr))
    return {'ierr': ierr.value, 'init_FT': float(new[8]), 'init_MT': float(new[10]),
            'delta0_FT_mm': float(new[12]), 'delta0_MT_mm': float(new[14]),
            'sigma0_FT_mpa': float(new[16]), 'sigma0_MT_mpa': float(new[18]),
            'kappa_MT_mm': float(new[6]), 'tau_mpa': float(s[shear_index])}


def run(core, run_dir):
    p = make_props("B", BASE_B, 0)
    lc = 2.23
    paths = {}
    # 1) 单增量反例（候选版：起始沿路径定位在 s*=gamma0/gamma，
    #    delta0=LC*gamma0、sigma0=S12=110；因 LC=2.23 仍固有不相容而报1033。
    #    修正前基于端点的文档反例（110.05/0.034565）存档于
    #    validation/reproduction_prefix_G1.json）。
    state, s, err = single_increment(core, p, 0.0155, lc)
    ply = strict_ply_forensics(core, 0.0155, lc, 3)
    sigma0 = ply['sigma0_MT_mpa']
    delta0 = ply['delta0_MT_mm']
    deltaf = 2*float(p[18])/sigma0 if sigma0 > 0 else None  # p[18]=PROPS(19)=GMT
    expect_d0 = lc*GAMMA0_12
    paths['single_increment_gamma12_0.0155'] = {
        'ierr': err, 'expected_ierr': 1033,
        'ply_forensics': ply,
        'sigma0_mt_mpa': sigma0, 'delta0_mt_mm': delta0,
        'deltaf_mt_mm': deltaf,
        'onset_located_at_path_point': bool(
            err == 1033 and ply['ierr'] == 33
            and abs(sigma0-S12) <= 1e-8
            and abs(delta0-expect_d0) <= 1e-12
            and deltaf is not None and abs(deltaf-2*0.5/S12) <= 1e-12),
    }
    # 2) 增量斜坡：逐接受点FI/起始量/错误。
    #    PROPS(10:16)=XT,XC,YT,YC,S12,S13,S23=[2080,1250,60,290,110,110,110]，
    #    故纯剪 gamma12 与 gamma23 的基体/纤维拉伸起始同为 S=110 -> gamma0=110/7100。
    #    0.006 与 0.01 为未起始对照（应无错误），0.0155 跨过 gamma0（应1033）。
    for gamma in (0.006, 0.01, 0.0155):
        increments = 20
        rows, error = ramp(core, p, gamma, lc, increments)
        key = f'ramp_gamma12_{gamma}'
        paths[key] = {
            'increments_requested': increments,
            'accepted_points': len(rows), 'ierr_final': error,
            'rows': rows,
        }
        with (run_dir/f'{key}.csv').open('w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
    # 2b) 细斜坡定位起始增量：首个端点 gamma>=gamma0 的增量应报1033。
    fine_rows, fine_error = ramp(core, p, 0.016, lc, 2000)
    crossing = next((r for r in fine_rows if r['gamma'] >= GAMMA0_12), None)
    fine_ok = bool(
        fine_error == 1033 and crossing is not None
        and fine_rows[crossing['increment']-1]['ierr'] == 1033
        and all(r['ierr'] == 0 for r in fine_rows[:crossing['increment']-1]))
    paths['fine_ramp_gamma12_onset_location'] = {
        'increments': 2000, 'ierr_final': fine_error,
        'first_crossing_increment': None if crossing is None else crossing['increment'],
        'gamma_at_error': None if crossing is None else crossing['gamma'],
        'gamma0': GAMMA0_12, 'onset_located_at_threshold': fine_ok,
    }
    # 3) 旧版静默跳过对照。
    paths['legacy_silent_skip_gamma12_0.0155'] = legacy_control(core, 0.0155, lc)
    # 4) 端点起始步长依赖：同终点 gamma=0.02、Lc=0.3（有效域内），1步 vs 20步。
    step_dependence = {}
    state, s, err = single_increment(core, p, 0.02, 0.3)
    step_dependence['1_step'] = {'delta0_mt_mm': float(state[14]),
                                 'sigma0_mt_mpa': float(state[18]),
                                 'd_mt': float(state[2]), 'ierr': err}
    rows, error = ramp(core, p, 0.02, 0.3, 20, shear=5)
    last = rows[-1]
    step_dependence['20_step'] = {'delta0_mt_mm': last['delta0_MT'],
                                  'sigma0_mt_mpa': last['sigma0_MT'],
                                  'd_mt': last['d_MT'], 'ierr': error}
    d1 = step_dependence['1_step']
    d20 = step_dependence['20_step']
    step_dependence['onset_differs'] = bool(
        abs(d1['delta0_mt_mm']-d20['delta0_mt_mm']) > 1e-9
        or abs(d1['sigma0_mt_mpa']-d20['sigma0_mt_mpa']) > 1e-6)
    paths['endpoint_onset_step_dependence'] = step_dependence
    # 5) 解析基准数值复核。
    analytic = {
        'gamma0_23shear': GAMMA0_23,
        'gamma0_reference': 110./7100.,
        'lc_crit_mt_pure_shear_mm': LC_CRIT_MT_SHEAR,
        'lc_crit_reference': 2*7100*0.5/110.**2,
        'note': 'Lc_crit 是纯剪基体拉伸路径的条件，不是全瓶网格建议。',
    }
    strict_1033_reproduced = bool(
        paths['single_increment_gamma12_0.0155']['onset_located_at_path_point']
        and paths['ramp_gamma12_0.0155']['ierr_final'] == 1033
        and paths['ramp_gamma12_0.006']['ierr_final'] is None
        and paths['ramp_gamma12_0.01']['ierr_final'] is None
        and paths['fine_ramp_gamma12_onset_location']['onset_located_at_threshold'])
    legacy_ok = paths['legacy_silent_skip_gamma12_0.0155']['silent_skip']
    evidence = {
        'stage': '1_reproduce', 'gate': 'G1',
        'build_note': '候选版重跑：起始量应沿路径定位（步长无关）；修正前的'
                      '端点起始证据存档于 validation/reproduction_prefix_G1.json。',
        'compiler': 'gfortran（版本见validation/environment.json与编译日志）',
        'dll_note': '本套件在唯一运行目录重新编译候选src/wcm_core.for与'
                    'baseline/hashin.for，未加载任何历史DLL。',
        'analytic': analytic,
        'paths': {k: ({kk: vv for kk, vv in v.items() if kk != 'rows'}
                      | {'rows': f'运行目录 {run_dir.name}/{k}.csv'})
                   if 'rows' in v else v
                   for k, v in paths.items()},
        'checks': {
            'strict_1033_reproduced_with_path_onset': strict_1033_reproduced,
            'onset_located_at_analytic_gamma0': paths[
                'fine_ramp_gamma12_onset_location']['onset_located_at_threshold'],
            'legacy_silent_skip_confirmed': legacy_ok,
            'onset_step_independent': not step_dependence['onset_differs'],
            'analytic_baselines_exact': bool(
                GAMMA0_23 == 110./7100. and LC_CRIT_MT_SHEAR == 7100./12100.),
        },
        'unresolved_1022': '见validation/forensics_1022.json：实际失败试算状态'
                           '不可从现有ODB/日志取得，保持UNRESOLVED_1022。',
    }
    evidence['status'] = 'PASS' if all(evidence['checks'].values()) else 'FAIL'
    return evidence
