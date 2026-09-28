# -*- coding: utf-8 -*-
"""WCM双族接口字典与独立的三维张量弹性参考。

A保留旧PROPS(1:22)，附加theta,wplus,update。
B保留旧PROPS(1:25)，附加theta,wplus,capf,capm,update。
前九项必须是原始单层工程常数，不是WCM Bin的等效常数。
强度、损伤参数仍保持原单层含义；无需按角度再乘cos或按族权重折减。
角度来自Beta/显式映射，不是Bin编号，不是AddRot。

update=0：保持原源码的分步更新。这不是当前增量全隐式损伤模型。
update=1：同一调用中反馈损伤；B使用中心差分切线，需UNSYMM。
B候选版起始量改为沿增量应变路径定位首次FI>=1的有效点（最早s0，
DEL/SIG非正者跳过），不再取增量端点；端点值仅驱动起始后的演化。
能量参数不合法时主动报错，不能靠改变G来掩盖网格过粗。

例：python tools/convert_wcm_inp.py --source ../Job-hashinnewconstant.inp --variant A
例：python tools/convert_wcm_inp.py --source ../Job-hashinnewcdm.inp --variant B
例：python tools/run_tests.py
例：python tools/abaqus_verify.py --suite elastic

旧项目损伤公式的物理假设、残余承载与能量非客观性不因角度修复消失。
此处不宣称复现Lin或Zhang完整模型，也不覆盖层间脱层和任意大应变。
"""
import numpy as np

SCHEMA = 202609
NSTATE = 84  # 核心计算的临时布局；不是新版 INP 的 *Depvar 数量。
MODE = {'A': 1, 'B': 2}
ERRORS = {
    1: '模型或材料参数数量错误（A25/B30）',
    2: '材料参数含NaN/Inf或过大值',
    3: '状态变量含NaN/Inf或过大值',
    4: '已废弃的材料版本标记错误',
    5: '旧SDV布局/重启动材料角度或版本不一致',
    6: '弹性模量、剪切模量或强度非正',
    7: 'theta或wplus超范围',
    8: 'update只能为0或1',
    9: '损伤参数、黏性或上限非法',
    10: '时间或特征长度非法',
    11: '应变含非有限数',
    12: '原始单层柔度不正定',
    13: '损伤历史超范围或起始标志不是0/1',
    21: 'FT起始等效量非正', 22: 'FC起始等效量非正',
    23: 'MT起始等效量非正', 24: 'MC起始等效量非正',
    31: 'FT的deltaf<=delta0', 32: 'FC的deltaf<=delta0',
    33: 'MT的deltaf<=delta0', 34: 'MC的deltaf<=delta0',
    41: 'FT未起始但增量起点已越限：缺少可定位起始历史',
    42: 'FC未起始但增量起点已越限：缺少可定位起始历史',
    43: 'MT未起始但增量起点已越限：缺少可定位起始历史',
    44: 'MC未起始但增量起点已越限：缺少可定位起始历史',
}


def error_text(code):
    family, reason = divmod(code, 1000)
    return f'族{family}: ' + ERRORS.get(reason, f'未知错误{reason}')


def elastic(constants):
    """从柔度求逆，与Fortran显式刚度表达式独立交叉验证。"""
    e1, e2, e3, v12, v13, v23, g12, g13, g23 = constants
    if min(e1, e2, e3, g12, g13, g23) <= 0:
        raise ValueError('单层模量必须为正')
    compliance = np.diag([1/e1, 1/e2, 1/e3, 1/g12, 1/g13, 1/g23])
    compliance[0, 1] = compliance[1, 0] = -v12/e1
    compliance[0, 2] = compliance[2, 0] = -v13/e1
    compliance[1, 2] = compliance[2, 1] = -v23/e2
    np.linalg.cholesky(compliance)
    return np.linalg.inv(compliance)


def tensor(v, strain=False):
    t = np.diag(v[:3])
    factor = 0.5 if strain else 1.0
    for value, (i, j) in zip(v[3:], ((0, 1), (0, 2), (1, 2))):
        t[i, j] = t[j, i] = value*factor
    return t


def vector(t, strain=False):
    factor = 2.0 if strain else 1.0
    return np.array([t[0, 0], t[1, 1], t[2, 2],
                     factor*t[0, 1], factor*t[0, 2], factor*t[1, 2]])


def rotation(theta):
    a = np.deg2rad(theta)
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])


def rotate_strain(v, theta):
    r = rotation(theta)
    return vector(r @ tensor(v, strain=True) @ r.T, strain=True)


def rotate_stress_back(v, theta):
    r = rotation(theta)
    return vector(r.T @ tensor(v) @ r)


def smeared_stiffness(constants, theta, weight=0.5):
    c = elastic(constants)
    result = np.zeros((6, 6))
    for j in range(6):
        e = np.eye(6)[j]
        for angle, w in ((theta, weight), (-theta, 1-weight)):
            result[:, j] += w*rotate_stress_back(c @ rotate_strain(e, angle), angle)
    return result


def make_props(variant, base, angle, weight=0.5, update=0,
               capf=1.0, capm=1.0):
    required = 22 if variant == 'A' else 25
    if len(base) != required:
        raise ValueError(f'{variant}原参数需{required}项，不接受静默补项')
    extra = [angle, weight]
    if variant == 'B':
        extra += [capf, capm]
    p = np.array(list(base)+extra+[update], dtype=float)
    validate_props(variant, p)
    return p


def validate_props(variant, p):
    if len(p) != (25 if variant == 'A' else 30) or not np.all(np.isfinite(p)):
        raise ValueError('材料参数数目错误或含非有限数')
    elastic(p[:9])
    if min(p[9:16]) <= 0:
        raise ValueError('强度必须为正')
    it = 22 if variant == 'A' else 25
    if not 0 <= p[it] <= 90 or not 0 <= p[it+1] <= 1:
        raise ValueError('theta/weight超范围')
    if p[-1] not in (0, 1):
        raise ValueError('update必须为0或1')
    if variant == 'A':
        if min(p[16:22]) < 0 or max(p[16:22]) > 1:
            raise ValueError('A损伤量/剪切权重必须在0到1内')
    else:
        if min(p[16:20]) <= 0 or min(p[21:25]) < 0 or max(p[21:23]) > 1:
            raise ValueError('B断裂能/黏性/剪切权重非法')
        if min(p[27:29]) <= 0 or max(p[27:29]) > 1:
            raise ValueError('B刚度折减上限必须大于0且不超过1')


def state_count(variant):
    return {'A': 20, 'B': 62}[variant]


def output_indices(variant):
    return tuple(range(1, {'A': 8, 'B': 10}[variant]+1))


def sdv_labels(variant):
    modes = ('FT', 'FC', 'MT', 'MC')
    labels = ['MAX_FI_'+m for m in modes]+['MAX_COMPLETE_F', 'MAX_COMPLETE_M']
    if variant == 'A':
        labels += ['MAX_USED_DF', 'MAX_USED_DM']
        local = ['INIT_'+m for m in modes]
    elif variant == 'B':
        labels += ['MAX_DVF', 'MAX_DVM', 'MAX_USED_DF', 'MAX_USED_DM']
        local = [prefix+'_'+m for prefix in ('D','KAPPA','INIT','DELTA0','SIGMA0','DV') for m in modes]
    else:
        raise ValueError(variant)
    labels += [family+'_'+label for family in ('PLUS','MINUS') for label in local]
    labels += ['THETA_DEG','WEIGHT_PLUS','MODEL_ID','SCHEMA']
    assert len(labels) == state_count(variant)
    return labels
