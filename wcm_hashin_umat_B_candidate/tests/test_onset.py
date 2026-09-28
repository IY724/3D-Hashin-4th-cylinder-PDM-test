# -*- coding: utf-8 -*-
"""阶段2起始定位验证套件（规划6.2/6.3/7，G2）。

独立参考：本模块按本构账本重新推导 FI/等效量公式并在应变路径上做
密集采样+二分求根，不调用被测 Fortran 的求根算法作答案。
覆盖：
- 单调/分支进入/多模式/非比例路径的起始量对比（相对容差按规划7）；
- 起始步长无关（1/20/200步同终点）；
- 错误语义：40+I 缺少历史、1023 无有效起始量仍拒绝；旧1022泊松反例
  改验应力等效FC起始，不冒充现场根因；基体不相容局部完成、纤维仍拒绝；
- 保留 deltaf=2G/sigma0 总断裂能口径，局部完成不宣称能量精确；
- 状态回滚：重复调用、拒绝试算、切线扰动、IADV=0 不污染历史；
- 先越限后回落：参考实现用于构造算例，再验证 Fortran 一致。
"""
import json
import sys
from pathlib import Path

import numpy as np

TOOLS = Path(__file__).resolve().parent.parent/'tools'
sys.path.insert(0, str(TOOLS))
from model import elastic, make_props  # noqa: E402

PLY = [141000., 11400., 11400., .28, .28, .28, 7100., 7100., 7100.]
STRENGTH = [2080., 1250., 60., 290., 110., 110., 110.]
BASE_B = PLY+STRENGTH+[133., 10., .5, 1.6, 1., .9, .5, 0., 0.]
C0 = elastic(PLY)
LC = 0.3
GAMMA0 = 110./7100.


def fi_ref(s, p):
    """独立重推的分支感知三维Hashin FI（与WCM_HASHIN同口径）。"""
    xt, xc, yt, yc, s12, s23 = p[9], p[10], p[11], p[12], p[13], p[15]
    fi = [0.0, 0.0, 0.0, 0.0]
    if s[0] >= 0:
        fi[0] = (s[0]/xt)**2+(s[3]/s12)**2+(s[4]/s12)**2
    else:
        fi[1] = (-s[0]/xc)**2
    q = s[1]+s[2]
    sh = (s[5]**2-s[1]*s[2])/s23**2+(s[3]/s12)**2+(s[4]/s12)**2
    if q >= 0:
        fi[2] = (q/yt)**2+sh
    else:
        fi[3] = ((yc/(2*s23))**2-1)*q/yc+q**2/(2*s23)**2+sh
    return fi


def eq_ref(e, s, lc, p=None):
    """独立等效量参考；FC用有效压应力/E11，不用总轴向压应变。"""
    e11 = PLY[0] if p is None else p[0]
    sh = e[3]**2+e[4]**2+e[5]**2
    shw = s[3]*e[3]+s[4]*e[4]+s[5]*e[5]
    d = [lc*np.sqrt(max(e[0], 0)**2+e[3]**2+e[4]**2),
         lc*max(-s[0], 0)/e11,
         lc*np.sqrt(max(e[1], 0)**2+max(e[2], 0)**2+sh),
         lc*np.sqrt(max(-e[1], 0)**2+max(-e[2], 0)**2+sh)]
    g = [0., 0., 0., 0.]
    if d[0] > 0:
        g[0] = lc*(max(s[0], 0)*max(e[0], 0)+s[3]*e[3]+s[4]*e[4])/d[0]
    if d[1] > 0:
        g[1] = max(-s[0], 0)
    if d[2] > 0:
        g[2] = lc*(max(s[1], 0)*max(e[1], 0)+max(s[2], 0)*max(e[2], 0)+shw)/d[2]
    if d[3] > 0:
        g[3] = lc*(max(-s[1], 0)*max(-e[1], 0)
                   +max(-s[2], 0)*max(-e[2], 0)+shw)/d[3]
    return d, g


def onset_ref(p, e0, e1, mode, samples=20001, lc=None):
    """独立路径搜索：最早满足 FI>=1 且 DEL>0、SIG>0 的 s；
    密集采样定位首个跨越区间后二分80次。返回 (s, delta0, sigma0) 或 None。"""
    lc = LC if lc is None else lc

    def pred(s):
        e = e0+s*(e1-e0)
        se = C0 @ e
        if fi_ref(se, p)[mode-1] < 1.0:
            return False
        d, g = eq_ref(e, se, lc, p)
        return d[mode-1] > 0 and g[mode-1] > 0
    grid = np.linspace(0., 1., samples)
    flags = [pred(s) for s in grid]
    hi = next((k for k, f in enumerate(flags) if f), None)
    if hi == 0:
        s = 0.0
    elif hi is None:
        return None
    else:
        lo, hi = grid[hi-1], grid[hi]
        for _ in range(80):
            mid = 0.5*(lo+hi)
            if pred(mid):
                hi = mid
            else:
                lo = mid
        s = hi
    e = e0+s*(e1-e0)
    se = C0 @ e
    d, g = eq_ref(e, se, lc, p)
    return s, d[mode-1], g[mode-1]


def fortran_onset(core, p, e0, e1, mode, lc=None):
    """单增量调用，返回族1的起始量与错误码。"""
    state, s, c, err = core.update('B', p, np.zeros(84), e0, e1-e0,
                                   length=LC if lc is None else lc,
                                   allow_error=True)
    return {'ierr': err,
            'init': float(state[7+mode]),          # NEW(8+mode)
            'delta0': float(state[11+mode]),       # NEW(12+mode)
            'sigma0': float(state[15+mode]),       # NEW(16+mode)
            'kappa': float(state[3+mode]),         # NEW(4+mode)
            'd': float(state[mode-1])}


def rel(a, b):
    return abs(a-b)/max(abs(b), 1e-30)


def mode_slots(mode):
    """族内指定模式的 d、kappa、init、delta0、sigma0、dv 六个历史槽。"""
    return np.arange(mode-1, 24, 4)


def inactive_histories_zero(state, active_modes):
    """两族未激活模式的全部历史必须为零，不能仅检查损伤值。"""
    return bool(all(np.all(state[ib+mode_slots(mode)] == 0.)
                    for ib in (0, 40) for mode in (1, 2, 3, 4)
                    if mode not in active_modes))


def update_with_old_probe(core, p, old, e0, de, lc):
    """绕过会复制OLD的Python包装，直接检查Fortran是否修改输入数组。"""
    import ctypes as ct
    props = np.array(p, dtype=np.float64, order='F', copy=True)
    e0, de = [np.array(a, dtype=np.float64, order='F', copy=True)
              for a in (e0, de)]
    new = old.copy()
    s, c = np.zeros(6), np.zeros((6, 6), order='F')
    mode, nprops, advance, err = map(ct.c_int, (2, len(props), 1, 0))
    dt, length = ct.c_double(.01), ct.c_double(lc)

    def ptr(a):
        return a.ctypes.data_as(ct.POINTER(ct.c_double))

    core.library.wcm_update_(ct.byref(mode), ptr(props), ct.byref(nprops),
                             ptr(old), ptr(e0), ptr(de), ct.byref(dt),
                             ct.byref(length), ct.byref(advance), ptr(new),
                             ptr(s), ptr(c), ct.byref(err))
    return new, s, c, err.value


def check_path(core, p, name, e0, e1, modes, tol=1e-8, cases=None):
    """单条路径：对指定模式对比 Fortran 与独立参考的起始量。"""
    row = {'name': name, 'e0': list(e0), 'e1': list(e1), 'modes': {}}
    ok = True
    for mode in modes:
        got = fortran_onset(core, p, e0, e1, mode)
        ref = onset_ref(p, e0, e1, mode)
        entry = {'fortran': got, 'reference': None if ref is None else
                 {'s': ref[0], 'delta0': ref[1], 'sigma0': ref[2]}}
        if ref is None:
            entry['pass'] = bool(got['ierr'] == 0 and got['init'] == 0)
        else:
            entry['pass'] = bool(
                got['ierr'] == 0 and got['init'] == 1
                and rel(got['delta0'], ref[1]) <= tol
                and rel(got['sigma0'], ref[2]) <= tol)
            entry['rel_delta0'] = rel(got['delta0'], ref[1])
            entry['rel_sigma0'] = rel(got['sigma0'], ref[2])
        row['modes'][mode] = entry
        ok = ok and entry['pass']
    row['pass'] = bool(ok)
    if cases is not None:
        cases.append(row)
    return row


def run(core, run_dir):
    p = make_props('B', BASE_B, 0)
    checks = {}
    cases = []

    # 1) 单调起始：拉/压用自由泊松方向（规划7：单轴应力必须用compliance
    #    构造，不得固定横向应变）；剪/混合用直接应变路径。
    z = np.zeros(6)
    dir1 = np.linalg.solve(C0, np.array([1., 0, 0, 0, 0, 0]))
    check_path(core, p, 'e1_ft_free_poisson', z, 2200*dir1, [1],
               cases=cases)
    check_path(core, p, 'e1_fc_free_poisson', z, -1300*dir1, [2],
               cases=cases)
    check_path(core, p, 'gamma23_mt', z, np.array([0, 0, 0, 0, 0, .02]), [3],
               cases=cases)
    check_path(core, p, 'gamma12_ft_mt', z, np.array([0, 0, 0, .02, 0, 0]), [1, 3],
               cases=cases)
    check_path(core, p, 'e2_mt_poisson', z, np.array([0, .02, 0, 0, 0, 0]), [3],
               cases=cases)
    # MC：等双轴横向压时 -S2*S3 项与 Q^2/(2S23)^2 相消，需叠加 23 剪切；
    # 该路径下 FT/FC（S1=0）与 MT（Q<0）分支均不激活，属 MC 单模式。
    check_path(core, p, 'mc_biaxial_shear', z,
               np.array([0, -.03, -.01, 0, 0, .03]), [4], cases=cases)
    check_path(core, p, 'mixed_e1_g12', z,
               np.array([.015, 0, 0, .01, 0, 0]), [1, 3], cases=cases)
    checks['monotonic_onsets_match_reference'] = all(c['pass'] for c in cases)

    # 1b) 固定横向应变时MT判据越限但DEL3=0：仍须报1023。
    #     这是无有效起始量，不属于已取得正起始量后的局部完成分支。
    _, _, _, err = core.update('B', p, np.zeros(84), z,
                               np.array([.02, 0, 0, 0, 0, 0]),
                               length=LC, allow_error=True)
    checks['fixed_lateral_uniaxial_mt_unseedable_1023'] = bool(err == 1023)

    # 2) 分支进入：S1 由负转正时剪切已越限（FT 分支入口）。
    branch = check_path(core, p, 'ft_branch_entry',
                        np.array([-.005, 0, 0, .005, 0, 0]),
                        np.array([.002, 0, 0, .02, 0, 0]), [1, 2, 3],
                        cases=cases)
    checks['branch_entry_handled'] = branch['pass']

    # 3) 解析锚点：gamma23 起始必须精确落在 gamma0=110/7100。
    got = fortran_onset(core, p, z, np.array([0, 0, 0, 0, 0, .02]), 3)
    checks['gamma23_analytic_anchor'] = bool(
        got['ierr'] == 0 and got['init'] == 1
        and abs(got['delta0']-LC*GAMMA0) <= 1e-12
        and abs(got['sigma0']-110.) <= 1e-8)

    # 4) 起始步长无关：1/20/200 步到同一终点。
    step_rows = {}
    d20 = None
    for n in (1, 20, 200):
        state = np.zeros(84)
        e0 = np.zeros(6)
        de = np.array([0, 0, 0, 0, 0, .02/n])
        for _ in range(n):
            state, s, c, err = core.update('B', p, state, e0, de,
                                           length=LC, allow_error=True)
            if err:
                break
            e0 = e0+de
        step_rows[n] = {'delta0': float(state[14]), 'sigma0': float(state[18]),
                        'd': float(state[2]), 'ierr': err}
        if n == 20:
            d20 = state.copy()
    checks['onset_step_independent'] = bool(
        all(r['ierr'] == 0 and abs(r['sigma0']-110.) <= 1e-8
            for r in step_rows.values())
        and abs(step_rows[1]['delta0']-LC*GAMMA0) <= 1e-12
        and abs(step_rows[20]['delta0']-LC*GAMMA0) <= 1e-12
        and abs(step_rows[200]['delta0']-LC*GAMMA0) <= 1e-12
        and abs(step_rows[1]['sigma0']-110.) <= 1e-8
        and abs(step_rows[1]['d']-step_rows[20]['d']) <= 1e-9
        and abs(step_rows[1]['d']-step_rows[200]['d']) <= 1e-9)

    # 5) 错误语义 40+I：未起始但 E0 已越限（缺少可定位历史）。
    #    MT：gamma23=0.02 使 FI3(E0) 越限；FC：S1=4512.8*(-0.28)=-1263.6。
    _, _, _, err = core.update('B', p, np.zeros(84),
                               np.array([0, 0, 0, 0, 0, .02]), np.zeros(6),
                               length=LC, allow_error=True)
    _, _, _, err2 = core.update('B', p, np.zeros(84),
                                np.array([0, -.28, 0, 0, 0, 0]), np.zeros(6),
                                length=LC, allow_error=True)
    checks['error_40_missing_history'] = bool(err == 1043 and err2 == 1042)

    # 6) 旧1022合成路径保留：轴向应变为零，泊松耦合仍使S11跨越-XC。
    #    起点MC早已越限，先从零加载建立真实历史，不能把1044混作FC错误。
    ep0 = np.array([0, -.24, 0, 0, 0, 0])
    ep1 = np.array([0, -.30, 0, 0, 0, 0])
    previous, _, _, prep_err = core.update('B', p, np.zeros(84), z, ep0,
                                          length=LC, allow_error=True)
    snapshot = previous.copy()
    current, _, _, err = core.update('B', p, previous, ep0, ep1-ep0,
                                     length=LC, allow_error=True)
    ref_fc = onset_ref(p, ep0, ep1, 2)
    d0_fc, s0_fc = LC*p[10]/p[0], p[10]
    df_fc = 2*p[17]/s0_fc
    kap_fc = LC*max(-(C0 @ ep1)[0], 0.)/p[0]
    expected_fc = df_fc*(kap_fc-d0_fc)/(kap_fc*(df_fc-d0_fc))
    legacy, _, _ = core.original('B', np.zeros(24), z, ep0, length=LC)
    legacy, _, _ = core.original('B', legacy, ep0, ep1-ep0, length=LC)
    legacy_unseedable = bool(legacy[9] == 1 and legacy[13] == 0.
                             and legacy[17] == 0. and legacy[1] == 0.)
    checks['poisson_fc_stress_equivalent_onset'] = bool(
        prep_err == 0 and err == 0 and np.array_equal(previous, snapshot)
        and ref_fc is not None and rel(ref_fc[1], d0_fc) <= 1e-8
        and rel(ref_fc[2], s0_fc) <= 1e-8 and df_fc > d0_fc
        and 0. < expected_fc < 1. and legacy_unseedable
        and inactive_histories_zero(current, (2, 4))
        and all(current[ib+9] == 1
                and rel(current[ib+13], d0_fc) <= 1e-8
                and rel(current[ib+17], s0_fc) <= 1e-8
                and rel(current[ib+5], kap_fc) <= 1e-8
                and abs(current[ib+1]-expected_fc) <= 1e-10
                and current[ib+21] == current[ib+1]
                and np.array_equal(current[ib+mode_slots(4)],
                                   snapshot[ib+mode_slots(4)])
                for ib in (0, 40)))
    poisson_evidence = {
        'prepare_ierr': prep_err, 'ierr': err, 'reference': ref_fc,
        'delta0_fc': float(current[13]), 'sigma0_fc': float(current[17]),
        'deltaf_fc': float(df_fc), 'expected_d_fc': float(expected_fc),
        'legacy_unseedable_confirmed': legacy_unseedable,
        'note': '旧应变等效量为零、旧版静默冻结；历史严格版对应1022机制，'
                '不据此认定现场根因已解决。'}

    # 7) 同一12剪切路径：小长度正常软化，大长度仅MT局部完成。
    #    FT也会起始，但必须按自己的有效软化区间演化，不能被MT置1污染。
    e12 = np.array([0, 0, 0, .0155, 0, 0])
    got = fortran_onset(core, p, z, e12, 3)
    lc_large = 2.23
    state, s, c, err = core.update('B', p, np.zeros(84), z, e12,
                                   length=lc_large, allow_error=True)
    d0, sigma0, kap = lc_large*GAMMA0, 110., lc_large*e12[3]
    df_mt, df_ft = 2*p[18]/sigma0, 2*p[16]/sigma0
    expected_ft = df_ft*(kap-d0)/(kap*(df_ft-d0))
    checks['matrix_local_completion_with_correct_onset'] = bool(
        err == 0 and got['ierr'] == 0 and got['init'] == 1
        and abs(got['delta0']-LC*GAMMA0) <= 1e-12
        and abs(got['sigma0']-sigma0) <= 1e-8
        and df_mt < d0 < df_ft and 0. < expected_ft < 1.
        and np.all(np.isfinite(s)) and np.all(np.isfinite(c))
        and inactive_histories_zero(state, (1, 3))
        and all(state[ib+2] == 1. and state[ib+22] == 1.
                and state[ib+10] == 1. and state[ib+8] == 1.
                and abs(state[ib+14]-d0) <= 1e-12
                and abs(state[ib+18]-sigma0) <= 1e-8
                and abs(state[ib+12]-d0) <= 1e-12
                and abs(state[ib+16]-sigma0) <= 1e-8
                and abs(state[ib+6]-kap) <= 1e-12
                and abs(state[ib]-expected_ft) <= 1e-10
                and state[ib+20] == state[ib]
                for ib in (0, 40)))

    # 8) 状态回滚（规划7回滚行）。
    st1, s1, c1, _ = core.update('B', p, np.zeros(84), z,
                                 np.array([0, 0, 0, 0, 0, .02]), length=LC)
    st2, s2, c2, _ = core.update('B', p, np.zeros(84), z,
                                 np.array([0, 0, 0, 0, 0, .02]), length=LC)
    pure = bool(np.array_equal(st1, st2) and np.array_equal(s1, s2)
                and np.array_equal(c1, c2))
    # FC自由泊松路径在LC=2.23下确有deltaf<=delta0，仍须拒绝1032。
    # 材料卡和长度保持不变，回退小步只到尚未起始的区域，不降低拒绝门槛。
    e0 = -100.*dir1
    accepted, _, _, prep_err = core.update('B', p, np.zeros(84), z, e0,
                                           length=2.23, allow_error=True)
    old_snapshot = accepted.copy()
    trial_old = accepted.copy()
    dirty, _, _, rejected_err = update_with_old_probe(
        core, p, trial_old, e0, -1300.*dir1-e0, 2.23)
    old_unchanged = bool(np.array_equal(trial_old, old_snapshot))
    rejected_unchanged = bool(np.array_equal(dirty, old_snapshot))
    fc_d0, fc_df = 2.23*p[10]/p[0], 2*p[17]/p[10]
    rollback = bool(prep_err == 0 and rejected_err == 1032
                    and fc_df < fc_d0 and old_unchanged and rejected_unchanged)
    clean = old_snapshot.copy()
    for target in (-200., -300., -400.):
        endpoint = target*dir1
        clean, sc, cc, ec = core.update('B', p, clean, e0, endpoint-e0,
                                        length=2.23, allow_error=True)
        dirty, sd, cd, ed = core.update('B', p, dirty, e0, endpoint-e0,
                                        length=2.23, allow_error=True)
        rollback = bool(rollback and ec == 0 and ed == 0
                        and np.array_equal(clean, dirty)
                        and np.array_equal(sc, sd) and np.array_equal(cc, cd))
        e0 = endpoint
    # IADV=0 不推进历史（状态需带模式戳，否则触发错误5）。
    st0 = np.zeros(84)
    st0[0] = .3
    st0[80:84] = [p[25], p[26], 2, 202609]
    out, _, _, _ = core.update('B', p, st0, z, np.ones(6)*1e-5, advance=0)
    noadv = bool(np.array_equal(out, st0))
    # UPDATE=1 切线不改变历史（应力槽随UPDATE定义本就不同，只比历史槽位）。
    p1 = make_props('B', BASE_B, 0, update=1)
    sa, _, _, _ = core.update('B', p1, np.zeros(84), z,
                              np.array([0, 0, 0, 0, 0, .02]), length=LC)
    p0 = make_props('B', BASE_B, 0)
    sb, _, _, _ = core.update('B', p0, np.zeros(84), z,
                              np.array([0, 0, 0, 0, 0, .02]), length=LC)
    hist = np.r_[0:24, 40:64, 80:84]
    tangent_clean = bool(np.array_equal(sa[hist], sb[hist]))
    checks['state_rollback'] = bool(pure and rollback and noadv
                                    and tangent_clean)
    checks['rollback_detail'] = {'pure_function': pure,
                                 'rejected_trial_no_pollution': rollback,
                                 'rejected_ierr': rejected_err,
                                 'rejected_delta0_fc': float(fc_d0),
                                 'rejected_deltaf_fc': float(fc_df),
                                 'fortran_old_unchanged': old_unchanged,
                                 'rejected_new_equals_old': rejected_unchanged,
                                 'iadv0_no_advance': noadv,
                                 'tangent_no_state_change': tangent_clean}

    # 9) 先越限后回落：用参考实现搜索构造算例（E0 预起始、终点 FI<1、
    #    路径中存在有效起始），再验证 Fortran 一致。
    rng = np.random.default_rng(7)
    found = None
    for _ in range(4000):
        e0 = rng.uniform(-.008, .008, 6)
        e1 = e0+rng.uniform(-.01, .01, 6)
        ok0 = all(fi_ref(C0 @ e0, p)[m] < 1.0 for m in range(4))
        if not ok0:
            continue
        fi_end = fi_ref(C0 @ e1, p)
        for mode in (1, 2, 3, 4):
            if fi_end[mode-1] >= 1.0:
                continue
            r = onset_ref(p, e0, e1, mode)
            if r is not None:
                found = (mode, e0.copy(), e1.copy(), r)
                break
        if found:
            break
    if found:
        mode, e0, e1, r = found
        got = fortran_onset(core, p, e0, e1, mode)
        overshoot = {'mode': mode, 'e0': list(e0), 'e1': list(e1),
                     'reference_s': r[0], 'reference_delta0': r[1],
                     'reference_sigma0': r[2],
                     'fortran_init': got['init'],
                     'fortran_delta0': got['delta0'],
                     'fortran_sigma0': got['sigma0'],
                     'ierr': got['ierr'],
                     'pass': bool(got['ierr'] == 0 and got['init'] == 1
                                  and rel(got['delta0'], r[1]) <= 1e-8
                                  and rel(got['sigma0'], r[2]) <= 1e-8)}
        checks['overshoot_then_return'] = overshoot['pass']
    else:
        overshoot = {'pass': False, 'note': '未找到先越限后回落算例'}
        checks['overshoot_then_return'] = False

    # 10) 0度双族对称：族1/族2起始量一致。
    sym = bool(np.array_equal(d20[0:24], d20[40:64]))
    checks['zero_degree_family_symmetry'] = sym

    with (run_dir/'onset_cases.json').open('w', encoding='utf-8') as f:
        json.dump(cases, f, ensure_ascii=False, indent=1)
    evidence = {
        'stage': '2_onset', 'gate': 'G2',
        'reference': '本模块 fi_ref/eq_ref/onset_ref 独立重推公式，'
                     '密集采样+二分，不调用被测Fortran算法。',
        'lc_mm': LC, 'gamma0': GAMMA0,
        'step_independence': step_rows,
        'overshoot_case': overshoot,
        'poisson_fc': poisson_evidence,
        'matrix_policy': 'deltaf=2G/sigma0不变；不相容MT局部完成，'
                         '非能量精确软化；非零黏性递推由matrix套件逐步核对。',
        'a_regression_note': 'A版回归由 tools/run_tests.py 的'
                             'test_original_formula_and_history_regression'
                             '（A vs 原版逐增量 <=2e-10）覆盖。',
        'cases': cases,
        'checks': checks,
    }
    core_checks = ['monotonic_onsets_match_reference', 'branch_entry_handled',
                   'gamma23_analytic_anchor', 'onset_step_independent',
                   'fixed_lateral_uniaxial_mt_unseedable_1023',
                   'error_40_missing_history', 'poisson_fc_stress_equivalent_onset',
                   'matrix_local_completion_with_correct_onset', 'state_rollback',
                   'overshoot_then_return', 'zero_degree_family_symmetry']
    evidence['status'] = 'PASS' if all(checks[k] for k in core_checks) else 'FAIL'
    return evidence
