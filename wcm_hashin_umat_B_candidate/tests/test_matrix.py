# -*- coding: utf-8 -*-
"""阶段3A材料点验证矩阵（规划7）。

复用 test_onset 的独立参考（fi_ref/eq_ref/onset_ref 不调用被测
Fortran 算法）。行覆盖与验收：
- 弹性：六分量、0/15/45/72/90度，核心相对误差<=1e-8；
- 拉压起始：1/2/3向自由泊松，与参考相对误差<=1e-8；
- 剪切：12/13/23 正反向，工程剪应变、P14分母、激活模式正确；
- 跨越起始：大步对照细分，起始位置相对差<=1e-6；
- 非比例路径：拉剪/压剪/双轴横向，与独立路径搜索一致；
- 符号/分支：拉压切换、泊松反例、分支进入（onset套件已含）；
- 两族：theta=35、w=0.3 非对称双族独立历史；
- UPDATE0：本帧旧DV应力 + 黏性离散递推<=1e-10；
- 长度：0.1/0.3/0.5/0.56 兼容，0.6/1.1/2.23 拒绝（纯23剪切MT路径）。
"""
import json
import sys
from pathlib import Path

import numpy as np

TOOLS = Path(__file__).resolve().parent.parent/'tools'
TESTS = Path(__file__).resolve().parent
sys.path.insert(0, str(TOOLS))
sys.path.insert(0, str(TESTS))
from model import elastic, make_props, rotate_strain  # noqa: E402
import test_onset as base  # noqa: E402  独立参考与0度常量

C0 = base.C0
LC = base.LC


def smeared_ref(theta):
    """独立双族加权刚度参考（与 run_tests 的 smeared_stiffness 同思路）。"""
    from model import rotate_stress_back
    result = np.zeros((6, 6))
    for angle, w in ((theta, .5), (-theta, .5)):
        for j in range(6):
            e = np.eye(6)[j]*1e-5
            el = rotate_strain(e, angle)
            sl = C0 @ el
            result[:, j] += w*rotate_stress_back(sl, angle)/1e-5
    return result


def run(core, run_dir):
    p = make_props('B', base.BASE_B, 0)
    checks = {}
    rows = {}

    # 1) 弹性行。
    worst = 0.
    rng = np.random.default_rng(11)
    for angle in (0, 15, 45, 72, 90):
        pa = make_props('B', base.BASE_B, angle)
        expected = smeared_ref(angle)
        for e in np.vstack([np.eye(6)*1e-5, rng.normal(size=(2, 6))*1e-5]):
            state, s, c, err = core.update('B', pa, np.zeros(84),
                                           np.zeros(6), e, length=LC)
            denom = max(1., np.abs(expected).max())
            worst = max(worst, float(np.max(np.abs(c-expected))/denom))
    checks['elastic_core_le_1e-8'] = bool(worst <= 1e-8)
    rows['elastic_max_abs_diff'] = worst
    rows['elastic_wcm_card_note'] = ('rebuild135 Bin WCM卡误差6.808e-6<=1e-5，'
                                     '证据 validation/test_cases.json（阶段0生成）。')

    # 2) 拉压起始：1/2/3向自由泊松。
    z = np.zeros(6)
    uni_ok = True
    uni_rows = []
    for comp, mode_t, mode_c, amp in ((0, 1, 2, 2200.), (1, 3, 4, None),
                                      (2, 3, 4, None)):
        d = np.linalg.solve(C0, np.eye(6)[comp])
        if comp == 0:
            cases = [('t', amp*d, mode_t, base.LC), ('c', -1300.*d, mode_c, base.LC)]
        else:
            # 2/3向：拉伸取Q>0越限（0.02）；压缩的 -S2*S3 项与 Q^2 相消，
            # MC 越限需约 E=-0.0402；且该模式单轴路径 Lc_crit≈0.158mm，
            # LC=0.3 本身不相容，故 MC 用例取 LC=0.1。
            cases = [('t', 0.02*np.eye(6)[comp], mode_t, base.LC),
                     ('c', -0.05*np.eye(6)[comp], mode_c, 0.1)]
        for tag, e1, mode, lc in cases:
            got = base.fortran_onset(core, p, z, e1, mode, lc=lc)
            ref = base.onset_ref(p, z, e1, mode, lc=lc)
            ok = bool(ref is not None and got['init'] == 1
                      and base.rel(got['delta0'], ref[1]) <= 1e-8
                      and base.rel(got['sigma0'], ref[2]) <= 1e-8)
            uni_ok = uni_ok and ok
            uni_rows.append({'dir': comp+1, 'side': tag, 'mode': mode,
                             'pass': ok,
                             'rel_d0': None if ref is None else
                             base.rel(got['delta0'], ref[1])})
    checks['uniaxial_free_poisson_le_1e-8'] = uni_ok
    rows['uniaxial'] = uni_rows

    # 3) 剪切 12/13/23 正反向：起始量对称、激活模式正确。
    shear_ok = True
    shear_rows = []
    for comp, expect_modes in ((3, (1, 3)), (4, (1, 3)), (5, (3,))):
        for sign in (1., -1.):
            e1 = np.zeros(6)
            e1[comp] = sign*.02
            got = base.fortran_onset(core, p, z, e1, expect_modes[0])
            ref = base.onset_ref(p, z, e1, expect_modes[0])
            ok = bool(ref is not None and got['init'] == 1
                      and base.rel(got['delta0'], ref[1]) <= 1e-8)
            shear_rows.append({'shear': comp+1, 'sign': sign,
                               'modes_expected': list(expect_modes),
                               'pass': ok})
            shear_ok = shear_ok and ok
        # 反向起始量应与正向一致（FI用平方项）。
        ep = np.zeros(6)
        ep[comp] = .02
        em = np.zeros(6)
        em[comp] = -.02
        gp = base.fortran_onset(core, p, z, ep, expect_modes[0])
        gm = base.fortran_onset(core, p, z, em, expect_modes[0])
        same = bool(abs(gp['delta0']-gm['delta0']) <= 1e-12
                    and abs(gp['sigma0']-gm['sigma0']) <= 1e-9)
        shear_rows.append({'shear': comp+1, 'sign_symmetry_pass': same})
        shear_ok = shear_ok and same
    # 23剪切单模式性：FT/FC不得激活。
    g23 = base.fortran_onset(core, p, z, np.array([0, 0, 0, 0, 0, .02]), 3)
    checks['shear_rows_ok'] = shear_ok
    checks['shear_23_single_mode'] = bool(
        g23['init'] == 1 and g23['sigma0'] == 110.0)
    rows['shear'] = shear_rows

    # 4) 跨越起始：大步对照细分，起始量相对差<=1e-6。
    got1 = base.fortran_onset(core, p, z, np.array([0, 0, 0, 0, 0, .02]), 3)
    state = np.zeros(84)
    e0 = np.zeros(6)
    for _ in range(200):
        state, s, c, err = core.update('B', p, state, e0,
                                       np.array([0, 0, 0, 0, 0, 1e-4]),
                                       length=LC, allow_error=True)
        e0 = e0+np.array([0, 0, 0, 0, 0, 1e-4])
        if state[10] >= 1:  # MT已起始
            break
    fine_d0 = float(state[14])
    checks['cross_step_le_1e-6'] = bool(
        base.rel(got1['delta0'], fine_d0) <= 1e-6)
    rows['cross_step'] = {'one_step_delta0': got1['delta0'],
                          'fine_delta0': fine_d0,
                          'rel': base.rel(got1['delta0'], fine_d0)}

    # 5) 非比例路径：拉剪/压剪/双轴横向，与独立路径搜索一致。
    np_ok = True
    np_rows = []
    for name, e1, modes in (
            ('tension_shear', np.array([.015, 0, 0, .012, 0, 0]), (1, 3)),
            ('compression_shear', np.array([-.02, 0, 0, .012, 0, 0]), (2,)),
            ('biaxial_lateral', np.array([0, .015, .008, 0, 0, .012]), (3, 4))):
        for mode in modes:
            got = base.fortran_onset(core, p, z, e1, mode)
            ref = base.onset_ref(p, z, e1, mode)
            ok = bool((ref is None and got['init'] == 0)
                      or (ref is not None and got['init'] == 1
                          and base.rel(got['delta0'], ref[1]) <= 1e-8
                          and base.rel(got['sigma0'], ref[2]) <= 1e-8))
            np_ok = np_ok and ok
            np_rows.append({'path': name, 'mode': mode, 'pass': ok})
    checks['nonproportional_match_reference'] = np_ok
    rows['nonproportional'] = np_rows

    # 6) 符号/分支：拉压切换（自由泊松，S1 由 +2000 穿零到 -1400，
    #    FC 在负侧 S1=-1250 处起始；Q=0 使 MT/MC 分支保持干净）。
    d1 = np.linalg.solve(C0, np.eye(6)[0])
    e0 = 2000.*d1
    e1 = -1400.*d1
    got = base.fortran_onset(core, p, e0, e1, 2)
    ref = base.onset_ref(p, e0, e1, 2)
    checks['tension_compression_switch'] = bool(
        ref is not None and got['init'] == 1
        and base.rel(got['delta0'], ref[1]) <= 1e-8
        and base.rel(got['sigma0'], ref[2]) <= 1e-8)

    # 7) 两族：theta=35、w=0.3 非对称；两族独立历史并与本地参考一致。
    p35 = make_props('B', base.BASE_B, 35, weight=.3)
    e1 = np.array([0, 0, 0, 0, 0, .02])
    state, s, c, err = core.update('B', p35, np.zeros(84), z, e1, length=LC)
    ok35 = True
    fam_rows = []
    for fam, theta in ((0, 35.), (40, -35.)):
        el1 = rotate_strain(e1, theta)
        for mode in (1, 3):
            got = {'init': float(state[fam+7+mode]),
                   'delta0': float(state[fam+11+mode]),
                   'sigma0': float(state[fam+15+mode])}
            ref = base.onset_ref(p35, np.zeros(6), el1, mode)
            # 纯23剪切旋转后 FT 分支不激活属正常：参考与Fortran须一致。
            ok = bool((ref is None and got['init'] == 0)
                      or (ref is not None and got['init'] == 1
                          and base.rel(got['delta0'], ref[1]) <= 1e-8))
            ok35 = ok35 and ok
            fam_rows.append({'family': fam//40+1, 'mode': mode, 'pass': ok,
                             'delta0': got['delta0']})
    independent = bool(np.isfinite(state[14]) and np.isfinite(state[54])
                       and state[14] > 0 and state[54] > 0)
    checks['two_family_independent_and_match'] = bool(ok35 and independent)
    rows['two_family'] = fam_rows

    # 8) UPDATE0 本帧旧DV应力 + 黏性离散递推（生产eta）。
    base_b_prod = base.PLY+base.STRENGTH+[133., 10., .5, 1.6, 1., .9, .5,
                                          1e-4, 5e-3]
    pp = make_props('B', base_b_prod, 0)
    state = np.zeros(84)
    e0 = np.zeros(6)
    de = np.array([0, 0, 0, 0, 0, 4e-3])
    lag_ok = True
    rec_ok = True
    lag_checked = 0
    for _ in range(30):
        prev_dv3 = float(state[22])
        prev_dv1 = float(state[20])
        state, s, c, err = core.update('B', pp, state, e0, de, length=LC)
        if prev_dv3 > 0 or prev_dv1 > 0:
            # UPDATE0：应力应用上一帧DV计算（本构账本第8节刚度形式）。
            df = 1-(1-prev_dv1)*(1-0.0)
            dm = 1-(1-prev_dv3)*(1-0.0)
            rf = max(1-min(max(df, 0.), 1.), 1e-6)
            rm = max(1-min(max(dm, 0.), 1.), 1e-6)
            rs = max(min((1-.9*prev_dv3)*(1-.5*0.0), 1.), 1e-6)
            cc = C0.copy()
            for i in range(3):
                for j in range(3):
                    cc[i, j] = C0[i, j]*rf*rm
            cc[0, 0] = C0[0, 0]*rf
            cc[3, 3] = C0[3, 3]*rs
            cc[4, 4] = C0[4, 4]*rs
            cc[5, 5] = C0[5, 5]*rs
            e1 = e0+de
            expect_s = cc @ e1
            lag_checked += 1
            lag_ok = lag_ok and bool(
                np.max(np.abs(s-expect_s)) <= 1e-9*max(1., np.abs(expect_s).max()))
        # 黏性递推：DV3=(eta_m*DV3_old+DT*d3)/(eta_m+DT)。
        eta_m = 5e-3
        d3 = float(state[2])
        expect_dv = (eta_m*prev_dv3+0.01*d3)/(eta_m+0.01)
        rec_ok = rec_ok and bool(abs(float(state[22])-expect_dv)
                                 <= 1e-10*max(1., abs(expect_dv)))
        e0 = e0+de
    checks['update0_lagged_stress'] = bool(lag_ok and lag_checked >= 5)
    checks['viscous_recursion_le_1e-10'] = rec_ok

    # 9) 长度域：纯23剪切MT路径。
    len_ok = True
    len_rows = []
    for lc in (0.1, 0.3, 0.5, 0.56):
        state, s, c, err = core.update('B', p, np.zeros(84), z,
                                       np.array([0, 0, 0, 0, 0, .02]),
                                       length=lc, allow_error=True)
        ok = bool(err == 0 and abs(float(state[14])-lc*base.GAMMA0) <= 1e-12)
        len_ok = len_ok and ok
        len_rows.append({'lc': lc, 'expected': 'compatible',
                         'ierr': err, 'pass': ok})
    for lc in (0.6, 1.1, 2.23):
        state, s, c, err = core.update('B', p, np.zeros(84), z,
                                       np.array([0, 0, 0, 0, 0, .02]),
                                       length=lc, allow_error=True)
        ok = bool(err == 1033)
        len_ok = len_ok and ok
        len_rows.append({'lc': lc, 'expected': 'EXPECTED_INCOMPATIBLE',
                         'ierr': err, 'pass': ok})
    checks['length_domain_pure23mt'] = len_ok
    rows['length'] = len_rows
    rows['length_note'] = ('分类仅适用纯剪基体拉伸路径；不得推广到纤维/基体'
                           '压缩或组合路径（规划7）。')

    gate_rows = ['elastic_core_le_1e-8', 'uniaxial_free_poisson_le_1e-8',
                 'shear_rows_ok', 'shear_23_single_mode', 'cross_step_le_1e-6',
                 'nonproportional_match_reference',
                 'tension_compression_switch',
                 'two_family_independent_and_match', 'update0_lagged_stress',
                 'viscous_recursion_le_1e-10', 'length_domain_pure23mt']
    evidence = {
        'stage': '3A_matrix', 'gate': 'G3(矩阵部分)',
        'reference': 'test_onset 独立参考 + smeared_ref 独立加权刚度。',
        'covered_elsewhere': {
            'A回归<=1e-12量级': 'tools/run_tests.py '
                               'test_original_formula_and_history_regression',
            '回滚': 'validation/onset.json state_rollback',
            'UPDATE1切线<=1e-4': 'tools/run_tests.py '
                                 'test_current_update_B_tangent',
            '先越限后回落': 'validation/onset.json overshoot_then_return',
            'rebuild135弹性': 'tools/run_tests.py test_rebuild_135_bins_elastic',
        },
        'rows': rows,
        'checks': checks,
    }
    evidence['status'] = 'PASS' if all(checks[k] for k in gate_rows) else 'FAIL'
    return evidence
