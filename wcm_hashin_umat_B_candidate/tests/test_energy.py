# -*- coding: utf-8 -*-
"""阶段3B能量闭合门（规划8，G3）。

解析负对照（规划8.2，待执行项在本套件落地）：
0度、纯23剪切、eta=0、UPDATE1、有效Lc、仅基体拉伸激活。当前B剪切
退化式 RS=(1-0.9*dmtv)*(1-0.5*dmcv)，仅 dmtv 激活时实际牵引
    tau = 0.1*未损伤弹性牵引 + 0.9*理想完全线性软化牵引
（恒等式，逐点核对）；到 deltaf 后剩余 10% 刚度的储能被扣除，
连续极限耗能 = 0.9*Gmt = 0.45 N/mm，而不是输入的 Gmt=0.5 N/mm。

账本（规划8.1）：每接受步 ΔW=0.5(σn+σn+1):Δε（工程剪应变功配对）；
Ψ=0.5 σ·ε（本构当前损伤下线性，等价 0.5 εᵀC(DV_used)ε）；
ΔD=ΔW-ΔΨ；裂带能 = Lc×体积耗散密度。0度双族相同，按 0.5/0.5 权重
相加后等于单族账本（逐帧核对族对称）。
放行（规划8.3）：细化后耗能相对误差<=2%、相邻细化差<=1%、应力曲线
归一化差<=1%、负耗能仅允许数值量级。
"""
import json
import sys
from pathlib import Path

import numpy as np

TOOLS = Path(__file__).resolve().parent.parent/'tools'
TESTS = Path(__file__).resolve().parent
sys.path.insert(0, str(TOOLS))
sys.path.insert(0, str(TESTS))
from model import make_props  # noqa: E402
import test_onset as base  # noqa: E402

G23 = 7100.
GMT = 0.5
SIGMA0 = 110.
LC = 0.3
GAMMA_END = 0.045            # > gamma_f = deltaf/Lc = 0.0303
DELTA0 = LC*base.GAMMA0
DELTAF = 2*GMT/SIGMA0
TARGET = 0.9*GMT             # 0.45 N/mm


def ramp_energy(core, base25, n, lc, update=1):
    """纯23剪切斜坡加载，返回逐步账本与摘要。"""
    p_u = make_props('B', base25, 0, update=update)
    state = np.zeros(84)
    e0 = np.zeros(6)
    de = np.array([0, 0, 0, 0, 0, GAMMA_END/n])
    W = 0.
    D = 0.
    min_step_D = 0.
    curve = []
    prev_tau = 0.
    prev_gam = 0.
    family_symmetric = True
    for i in range(1, n+1):
        state, s, c, err = core.update('B', p_u, state, e0, de,
                                       length=lc, allow_error=True)
        if err:
            break
        gam = GAMMA_END*i/n
        tau = float(s[5])
        dW = 0.5*(prev_tau+tau)*(gam-prev_gam)
        W += dW
        # ΔD 用本帧与上帧储能差：Ψ=0.5 σ·ε = 0.5 τ γ（纯23剪切）。
        psi_now = 0.5*tau*gam
        if i == 1:
            d_psi = psi_now
        else:
            d_psi = psi_now-curve[-1]['psi']
        dD = dW-d_psi
        D += dD
        min_step_D = min(min_step_D, dD)
        curve.append({'i': i, 'gamma': gam, 'tau': tau,
                      'd': float(state[2]), 'psi': psi_now, 'dW': dW,
                      'dD': dD})
        family_symmetric = family_symmetric and bool(
            np.array_equal(state[0:24], state[40:64]))
        prev_tau, prev_gam = tau, gam
        e0 = e0+de
    D_crack = lc*D
    return {'n': n, 'lc': lc, 'update': update, 'W_density': W,
            'psi_final': curve[-1]['psi'] if curve else None,
            'D_density': D, 'D_crack_N_per_mm': D_crack,
            'rel_error_vs_target': abs(D_crack-TARGET)/TARGET,
            'min_step_dD': min_step_D,
            'family_symmetric': family_symmetric,
            'curve': curve}


def traction_identity_check(curve):
    """逐点核对 tau = 0.1*G23*gamma + 0.9*tau_ideal(delta)。
    恒等式只在起始后（delta>=delta0）成立；未起始弹性段跳过。"""
    worst = 0.
    checked = 0
    for pt in curve:
        delta = LC*pt['gamma']
        if delta < DELTA0:
            continue
        tau_ideal = SIGMA0*(DELTAF-delta)/(DELTAF-DELTA0) if delta < DELTAF else 0.
        expect = 0.1*G23*pt['gamma']+0.9*max(tau_ideal, 0.)
        worst = max(worst, abs(pt['tau']-expect)/max(1., abs(expect)))
        checked += 1
    traction_identity_check.checked = checked
    return worst


def psi_explicit_check(curve):
    """Ψ = 0.5*G23*(1-0.9d)*γ² 显式重构对比（当前DV=d，eta=0）。"""
    worst = 0.
    for pt in curve:
        expect = 0.5*G23*(1-0.9*pt['d'])*pt['gamma']**2
        worst = max(worst, abs(pt['psi']-expect)/max(1e-12, abs(expect)))
    return worst


def run(core, run_dir):
    base25_eta0 = base.BASE_B             # 25项基参数，eta=0
    checks = {}
    rows = {}
    # 1) eta=0、UPDATE1、LC=0.3，三档细化。
    runs = {}
    for n in (100, 200, 400):
        runs[n] = ramp_energy(core, base25_eta0, n, LC, update=1)
        runs[n].pop('curve')
    errs = [runs[n]['rel_error_vs_target'] for n in (100, 200, 400)]
    refine_gap = max(abs(errs[0]-errs[1]), abs(errs[1]-errs[2]))
    full_curve = ramp_energy(core, base25_eta0, 400, LC, update=1)['curve']
    id_worst = traction_identity_check(full_curve)
    psi_worst = psi_explicit_check(full_curve)
    checks['update1_eta0_dissipation_le_2pct'] = bool(errs[-1] <= 0.02)
    checks['refinement_gap_le_1pct'] = bool(refine_gap <= 0.01)
    checks['traction_identity_le_1e-9'] = bool(id_worst <= 1e-9)
    rows['traction_checked_points'] = int(getattr(traction_identity_check, 'checked', 0))
    checks['psi_explicit_le_1e-10'] = bool(psi_worst <= 1e-10)
    checks['no_negative_dissipation'] = bool(
        runs[400]['min_step_dD'] >= -1e-6*max(runs[400]['W_density'], 1e-6))
    checks['family_symmetric_0deg'] = bool(runs[400]['family_symmetric'])
    rows['update1_eta0'] = runs
    rows['traction_identity_worst'] = id_worst
    rows['psi_explicit_worst'] = psi_worst
    rows['refinement_gap'] = refine_gap

    # 2) UPDATE0（eta=0）细步极限：滞后一帧，误差应随细化下降。
    u0 = {n: ramp_energy(core, base25_eta0, n, LC, update=0) for n in (100, 400)}
    for r in u0.values():
        r.pop('curve')
    checks['update0_converges_to_target'] = bool(
        u0[400]['rel_error_vs_target'] < u0[100]['rel_error_vs_target'])
    rows['update0_eta0'] = u0

    # 3) 生产黏性（eta_f=1e-4, eta_m=5e-3）UPDATE0：账本含黏性附加耗散，
    #    只要求无负耗散与递推闭合（递推已在matrix套件<=1e-10核对）。
    base25_prod = base.PLY+base.STRENGTH+[133., 10., .5, 1.6, 1., .9, .5,
                                          1e-4, 5e-3]
    prod = ramp_energy(core, base25_prod, 400, LC, update=0)
    prod.pop('curve')
    checks['production_eta_no_negative_dissipation'] = bool(
        prod['min_step_dD'] >= -1e-6*max(prod['W_density'], 1e-6))
    rows['production_eta_update0'] = prod

    # 4) 有效长度域上限 LC=0.56：解析极限仍应为 0.45。
    run56 = ramp_energy(core, base25_eta0, 200, 0.56, update=1)
    run56.pop('curve')
    checks['lc056_same_limit_045'] = bool(
        run56['rel_error_vs_target'] <= 0.02)
    rows['lc056_update1_eta0'] = run56

    # 5) 与输出SDV/SSE的关系说明（规划8.1）：SPD为有符号诊断，
    #    不能当单模式断裂能；生产SDV1-10不足以重建分模式能量。
    rows['notes'] = ('SPD/SSE为umat_entry有符号诊断输出，Abaqus级SPD核对属'
                     '阶段4 fixture（G4）；本套件用应力/应变完整历史记账。')

    gate = ['update1_eta0_dissipation_le_2pct', 'refinement_gap_le_1pct',
            'traction_identity_le_1e-9', 'psi_explicit_le_1e-10',
            'no_negative_dissipation', 'family_symmetric_0deg',
            'update0_converges_to_target',
            'production_eta_no_negative_dissipation', 'lc056_same_limit_045']
    evidence = {
        'stage': '3B_energy', 'gate': 'G3',
        'analytic_target_N_per_mm': TARGET,
        'delta0_mm': DELTA0, 'deltaf_mm': DELTAF,
        'gamma_f': DELTAF/LC,
        'interpretation': ('0.45=0.9*Gmt 是当前退化式(1-0.9dmt)决定的连续极限；'
                           'Gmt=0.5 的输入名义值与实际耗能差10%，属模型定义'
                           '层面（剪切通道共享），不是数值误差（规划8.2）。'),
        'rows': rows,
        'checks': checks,
    }
    evidence['status'] = 'PASS' if all(checks[g] for g in gate) else 'FAIL'
    with (run_dir/'energy_curve_400.csv').open('w', encoding='utf-8') as f:
        f.write('i,gamma,tau,d,psi,dW,dD\n')
        for pt in full_curve:
            f.write(f"{pt['i']},{pt['gamma']},{pt['tau']},{pt['d']},"
                    f"{pt['psi']},{pt['dW']},{pt['dD']}\n")
    return evidence
