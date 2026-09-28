# -*- coding: utf-8 -*-
"""Type-IV 氢气瓶建模插件 · 对话框设计样图 v4（两栏原创版式，基准=geometry_zhuning_joint_rebuild）。
渲染三种封头选项的联动效果各一张：
  hemispherical / semiellipsoidal / isotensoid
字段随选项切换：k 仅半椭球可用，α0 仅等张力可用，l2 在半球/等张力下自动计算（锁定）；
示意图封头曲线随选项更换（半球=圆弧、半椭球=提取的椭圆弧、等张力=等张力形示意曲线）。
轮廓取自 results/profile_extract.json（经体积校验的新模型真实子午线）。不需要 Abaqus。

运行：python scripts/draw_plugin_dialog_mock.py
输出：results/type4_plugin_mock_head_{hemispherical,semiellipsoidal,isotensoid}.png
     （半椭球版同时覆盖 results/type4_plugin_dialog_mock.png）
"""
import json
import math
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch, Circle

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

ROOT = Path(__file__).resolve().parent.parent
DATA = json.loads((ROOT / 'results' / 'profile_extract.json').read_text(encoding='utf-8'))

C = '#1a3fbf'
CW = '#263442'
UI = '#263442'
BOX = '#8b99a6'
GRAY = '#8a949e'

HEADS = ('hemispherical', 'semiellipsoidal', 'isotensoid')
HEAD_LABEL = {'hemispherical': 'Hemispherical', 'semiellipsoidal': 'Semi-ellipsoidal',
              'isotensoid': 'Isotensoid'}
HEAD_ICON = [lambda u: 1.35 * np.sqrt(np.maximum(0, 1 - u ** 2)),
             lambda u: 0.95 * np.sqrt(np.maximum(0, 1 - u ** 2)),
             lambda u: 1.15 * np.maximum(0, 1 - u ** 4) ** 0.38]


def dist(p, q):
    return math.hypot(p[0] - q[0], p[1] - q[1])


def chain_liner_segments():
    segs = DATA['liner']['segments']
    pool = list(segs)
    s0 = pool.pop(0)
    chain, start, end = [s0], s0['A'], s0['B']
    while dist(start, end) > 1e-6:
        for j, t in enumerate(pool):
            if dist(t['A'], end) < 1e-6:
                chain.append(t); end = t['B']; pool.pop(j); break
            if dist(t['B'], end) < 1e-6:
                f = dict(t); f['A'], f['B'] = t['B'], t['A']
                if f['kind'] == 'arc':
                    f['u1'], f['u2'] = t['u2'], t['u1']
                chain.append(f); end = f['B']; pool.pop(j); break
        else:
            raise RuntimeError('chain break at %s' % (end,))
    return chain


def sample_segment(s):
    if s['kind'] == 'line':
        return [tuple(s['A']), tuple(s['B'])]
    u1, u2, c, R = s['u1'], s['u2'], s['center'], s['radius']
    n = max(8, int(abs(u2 - u1) * R / 0.4))
    return [(c[0] + R * math.cos(u1 + (u2 - u1) * k / float(n)),
             c[1] + R * math.sin(u1 + (u2 - u1) * k / float(n))) for k in range(n + 1)]


R_OUT, R_SP_MM, H_MM = 140.0, 90.0, 119.0     # 外径/分界半径/极孔面-赤道轴向距（=内核默认）


def _iso_z_table():
    """等张力封头 ODE（与内核 _iso_trajectory 同款 RK4），返回 [(r, z)]，z 距赤道。"""
    alpha0 = math.radians(15.0)
    r0 = R_OUT * math.sin(alpha0)
    dz = 0.25

    def acc(r, rp):
        tan2 = r0 * r0 / max(r * r - r0 * r0, 1e-9)
        return -(1.0 + rp * rp) * (2.0 - tan2) / r

    r, rp, z = R_OUT, 0.0, 0.0
    pts = [(r, z)]
    for _ in range(20000):
        k1r, k1p = rp, acc(r, rp)
        k2r, k2p = rp + 0.5 * dz * k1p, acc(r + 0.5 * dz * k1r, rp + 0.5 * dz * k1p)
        k3r, k3p = rp + 0.5 * dz * k2p, acc(r + 0.5 * dz * k2r, rp + 0.5 * dz * k2p)
        k4r, k4p = rp + dz * k3p, acc(r + dz * k3r, rp + dz * k3p)
        r += dz / 6.0 * (k1r + 2 * k2r + 2 * k3r + k4r)
        rp += dz / 6.0 * (k1p + 2 * k2p + 2 * k3p + k4p)
        z += dz
        if r > R_OUT + 1e-6:
            break
        pts.append((r, z))
        if rp >= 0.0 and z > 1.0:
            break
    return pts


_ISO_PTS = _iso_z_table()


def head_z(head, r):
    """封头子午线 z(r)：距赤道平面的轴向深度（公式与内核 _head_pts 一致）。"""
    if head == 'hemispherical':
        return math.sqrt(max(R_OUT * R_OUT - r * r, 0.0))
    if head == 'isotensoid':
        for (ra, za), (rb, zb) in zip(_ISO_PTS[:-1], _ISO_PTS[1:]):
            if ra >= r >= rb:
                t = (ra - r) / max(ra - rb, 1e-12)
                return za + t * (zb - za)
        return _ISO_PTS[-1][1]
    b = R_OUT / 2.0                                  # semiellipsoidal（k=2）
    return b * math.sqrt(max(1.0 - (r / R_OUT) ** 2, 0.0))


def z_sp_of(head):
    """l2 = 分界半径 r_sp 处的封头深度（内核 derived['l2'] 同定义）。"""
    return head_z(head, R_SP_MM)


def r_at_z(head, z):
    """head_z 的反函数：给定距赤道深度 z，求封头子午线半径 r。"""
    if head == 'hemispherical':
        return math.sqrt(max(R_OUT * R_OUT - z * z, 0.0))
    if head == 'isotensoid':
        for (ra, za), (rb, zb) in zip(_ISO_PTS[:-1], _ISO_PTS[1:]):
            if za <= z <= zb or zb <= z <= za:
                t = (z - za) / max(zb - za, 1e-12)
                return ra + t * (rb - ra)
        return _ISO_PTS[-1][0]
    b = R_OUT / 2.0
    return R_OUT * math.sqrt(max(1.0 - (z / b) ** 2, 0.0))


def dome_curve(head, p0, p1):
    """封头外子午线：赤道 (140,641) → 分界 (90, 641+z_sp)，内核公式采样。"""
    n = 48
    pts = []
    for k in range(n + 1):
        r = R_OUT + (R_SP_MM - R_OUT) * k / float(n)
        pts.append((r, 641.0 + head_z(head, r)))
    return pts


LINER_CHAIN = chain_liner_segments()
DOME_SEG = next(s for s in LINER_CHAIN if s['A'][1] > 640 and s['B'][1] > 694 and s['A'][0] > 139)
BOSS_LOOP = max(DATA['boss']['loops'], key=lambda lp: min(q[1] for q in lp))


def liner_half_polyline(head):
    """+Y 半环采样折线（截断 y<430），封头段与分界壁按封头类型替换。"""
    pts = []
    z_sp = z_sp_of(head)
    for s in LINER_CHAIN:
        if s is DOME_SEG:
            sp = dome_curve(head, s['A'], s['B'])
        elif head != 'semiellipsoidal' and s['A'][0] == 90.0 and s['A'][1] > 690.0:
            sp = [(90.0, 641.0 + z_sp), tuple(s['B'])]   # 分界壁随 z_sp 伸长
        else:
            sp = sample_segment(s)
        pts += sp if not pts else sp[1:]
        if pts[-1][1] < 431.0:
            break
    out = [pts[0]]
    for q in pts[1:]:
        if dist(q, out[-1]) > 1e-9:
            out.append(q)
    return out


def render(head):
    fig, ax = plt.subplots(figsize=(11.6, 7.2))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 61)
    ax.set_aspect('equal')
    ax.axis('off')

    ax.add_patch(Rectangle((0.5, 0.5), 99, 60, fc='white', ec=BOX, lw=1.2))
    ax.add_patch(Rectangle((0.5, 57.6), 99, 2.9, fc='#dfe5ea', ec=BOX, lw=1.2))
    ax.text(2.2, 59.05, 'Type-IV Hydrogen Vessel  —  Liner & Boss Generator', fontsize=11,
            weight='bold', color=UI, va='center')

    # ---------- 左栏 ----------
    ax.text(3.0, 51.9, 'Inputs', fontsize=9.5, weight='bold', color=UI)
    ax.add_patch(Rectangle((2, 19.5), 38, 30.6, fc='#f4f6f8', ec=BOX, lw=1.0))

    def group(y, name):
        ax.text(3.2, y, name, fontsize=8.3, weight='bold', color=UI, va='center')
        ax.plot([3.2, 38.8], [y - 0.78] * 2, color='#c3ccd4', lw=0.8)

    def mfield(x, y, sym, val, lw=3.0, w=6.6, off=False):
        col = GRAY if off else UI
        ax.text(x, y, sym, fontsize=8.3, color=col, va='center')
        ax.add_patch(Rectangle((x + lw, y - 0.92), w, 1.95,
                               fc='#f0f0f0' if off else 'white', ec=BOX, lw=0.9))
        ax.text(x + lw + w - 0.4, y, val, fontsize=7.8, color=GRAY if off else '#555f6a',
                va='center', ha='right')

    group(49.2, 'Head type')
    for k, key in enumerate(HEADS):
        yc = 47.5 - k * 1.58
        sel = key == head
        ax.add_patch(Circle((3.7, yc + 0.5), 0.4, fc=C if sel else 'white', ec=BOX, lw=0.9))
        u = np.linspace(-1, 1, 30)
        ax.plot(4.9 + 1.5 * u, yc + 0.32 + HEAD_ICON[k](u), color=C if sel else GRAY, lw=1.2)
        ax.text(8.1, yc + 0.52, HEAD_LABEL[key], fontsize=8, color=C if sel else '#5a646e',
                va='center', weight='bold' if sel else 'normal')
    mfield(3.4, 42.4, '$k$', '2.0' if head == 'semiellipsoidal' else '—',
           off=head != 'semiellipsoidal')
    mfield(20.0, 42.4, '$\\alpha_0$', '55' if head == 'isotensoid' else '—',
           off=head != 'isotensoid')
    note = {'hemispherical': 'hemisphere: head depth $l_2$ computed from R, field locked',
            'semiellipsoidal': 'axis ratio k drives the dome; fields switch with head type',
            'isotensoid': 'Isotensoid: meridian by equal-tension ODE (uses $r_0$); $l_2$ auto'}[head]
    ax.text(3.4, 40.7, note, fontsize=6.8, color='#5a646e')

    l2_val = '%.1f' % z_sp_of(head)                  # 内核 derived['l2'] 同值
    l2_off = True                                    # l2 恒为自动计算

    group(38.9, 'Barrel  (cylindrical section)')
    for i, (s, v, lw) in enumerate([('$R$', '135', 2.4), ('$l_1$', '1282', 3.2), ('$t_1$', '5', 3.2)]):
        mfield(3.2 + i * 11.8, 37.2, s, v, lw=lw)

    group(35.3, 'Dome & liner shoulder')
    for i, (s, v, lw, off) in enumerate([('$l_2$', l2_val, 3.0, l2_off),
                                         ('$r_{so}$', '30', 4.8, False),
                                         ('$r_{si}$', '28', 4.8, False)]):
        mfield(3.2 + i * 11.8, 33.6, s, v, lw=lw, off=off)
    for i, (s, v, lw) in enumerate([('$r_c$', '50', 3.6), ('$l_{tip}$', '22', 4.8)]):
        mfield(3.2 + i * 16.9, 31.5, s, v, lw=lw, w=9.0)

    group(29.6, 'Polar joint  (liner-boss bond)')
    for i, (s, v, lw) in enumerate([('$r_{sp}$', '90', 4.8), ('$r_f$', '43', 3.6), ('$r_{if}$', '35', 4.8)]):
        mfield(3.2 + i * 11.8, 27.9, s, v, lw=lw)

    group(26.0, 'Polar boss  (metal insert)')
    for i, (s, v, lw) in enumerate([('$r_0$', '16', 3.0), ('$r_n$', '36', 3.6), ('$t_4$', '20', 3.2)]):
        mfield(3.2 + i * 11.8, 24.3, s, v, lw=lw)
    for i, (s, v, lw) in enumerate([('$r_{nf}$', '12', 4.8), ('$r_{bf}$', '6', 4.8)]):
        mfield(3.2 + i * 16.9, 22.2, s, v, lw=lw, w=9.0)

    # ---------- 右栏 ----------
    ax.text(43.5, 51.9, 'Schematic  (redrawn from inputs)', fontsize=9.5, weight='bold', color=UI)
    ax.add_patch(Rectangle((42, 19.5), 56, 30.6, fc='#ebebeb', ec=BOX, lw=1.0))

    X = lambda x: 43.5 + 0.163 * x     # x = 距极孔端面轴向距离 (= 760 - y)
    Y = lambda r: 22.5 + 0.150 * r
    endy = 760.0

    lp = liner_half_polyline(head)
    P = np.array([[X(endy - y), Y(r)] for r, y in lp])
    ax.add_patch(plt.Polygon(P, closed=True, fc='#cfe0f2', ec=CW, lw=1.5, zorder=2))
    P = np.array([[X(endy - y), Y(r)] for r, y in BOSS_LOOP])
    ax.add_patch(plt.Polygon(P, closed=True, fc='#f6dfc2', ec=CW, lw=1.5, zorder=2))

    ax.plot([42.6, 97.4], [Y(0)] * 2, color='#e08bb0', lw=1.0, ls=(0, (7, 3, 1, 3)), zorder=1)
    zx = X(760.0 - 430.0 + 10)
    zz = [(zx - .6, Y(127)), (zx + .6, Y(130.5)), (zx - .6, Y(133.5)),
          (zx + .6, Y(136.5)), (zx - .6, Y(139.5))]
    ax.plot([p[0] for p in zz], [p[1] for p in zz], color=CW, lw=1.0, zorder=3)

    DIM = dict(arrowstyle='<|-|>', color=C, lw=0.9, mutation_scale=7, shrinkA=0, shrinkB=0)

    def dim_h(x1, x2, y, text, ty=None):
        ax.annotate('', xy=(x1, y), xytext=(x2, y), arrowprops=DIM, zorder=5)
        ax.text((x1 + x2) / 2, ty if ty else y + 0.6, text, ha='center', va='bottom',
                color=C, fontsize=10, zorder=5)

    def dim_v(y1, y2, x, text):
        ax.annotate('', xy=(x, y1), xytext=(x, y2), arrowprops=DIM, zorder=5)
        ax.text(x + 0.6, (y1 + y2) / 2, text, ha='left', va='center', color=C,
                fontsize=10, rotation=90, zorder=5)

    def leader(xpt, ypt, xt, yt, text):
        ax.annotate(text, xy=(xpt, ypt), xytext=(xt, yt), fontsize=10, color=C,
                    ha='left', va='center', zorder=5,
                    arrowprops=dict(arrowstyle='->', color=C, lw=0.8, shrinkA=1, shrinkB=1))

    if head == 'semiellipsoidal':
        dim_h(X(H_MM - z_sp_of(head)), X(H_MM), 47.4, '$l_2$', ty=48.0)
    else:
        dim_h(X(H_MM - z_sp_of(head)), X(H_MM), 47.4, '$l_2$ (auto)', ty=48.0)
    dim_h(X(119), X(330), 49.1, '$l_1$', ty=49.7)
    dim_v(Y(0), Y(135), 95.0, '$R$')
    dim_h(X(0), X(330), 21.2, '$L$', ty=20.4)

    z_sp = z_sp_of(head)
    ax_sp = H_MM - z_sp                              # 分界平面轴向位置（距极孔端面）
    # 指引箭头：锚点取自 profile_extract.json 真实特征，路径两两不交叉
    leader(X(40), Y(16), 44.2, 23.8, '$r_0$')
    leader(X(30), Y(36), 45.0, 31.0, '$t_4$')
    leader(X(49.6), Y(39.5), 47.5, 33.5, '$r_{nf}$')
    leader(X((ax_sp + 82.0) / 2.0), Y(90), 48.5, 45.0, '$r_{sp}$')
    leader(X(111), Y(43.5), 67.5, 33.0, '$r_{if}$')
    leader(X(99), Y(r_at_z(head, H_MM - 99.0)), 61.5, 45.3, '$r_{so}$')
    leader(X(90), Y(71), 65.0, 38.5, '$r_{si}$')
    leader(X(129), Y(50), 71.0, 29.5, '$r_c$')
    leader(X(138), Y(43), 72.0, 27.0, '$r_f$')
    leader(X(141), Y(40), 75.0, 23.6, '$l_{tip}$')
    leader(X(116), Y(19), 70.0, 24.6, '$r_{bf}$')
    leader(X(250), Y(137.5), 88.0, 45.9, '$t_1$')

    # ---------- 底部 ----------
    ax.text(2.5, 16.6, 'Sector angle  $\\theta$:', fontsize=8.5, color=UI, va='center')
    ax.add_patch(Rectangle((13.2, 15.65), 7.0, 2.0, fc='white', ec=BOX, lw=0.9))
    ax.text(19.7, 16.65, '10', fontsize=8, color='#555f6a', va='center', ha='right')
    ax.text(20.9, 16.65, '°  (360 = full model)', fontsize=8, color='#5a646e', va='center')
    ax.text(46.0, 16.65, 'Reference standard:  GB/T 42612-2023', fontsize=8.5, color=UI, va='center')
    ax.text(2.5, 13.6, "Parts are created in the current model and named automatically:  'Liner' ,  'Boss'.\n"
                       "All dimensions in mm; values shown are an editable preset.",
            fontsize=7.5, color='#5a646e', va='center')

    for x0, w, name in [(68, 9, 'OK'), (78.5, 11, 'Apply'), (91, 8, 'Cancel')]:
        ax.add_patch(FancyBboxPatch((x0, 1.4), w, 2.9, boxstyle='round,pad=0.25',
                                    fc='#e4e9ee', ec=BOX, lw=1.0))
        ax.text(x0 + w / 2, 2.85, name, fontsize=9, color=UI, ha='center', va='center')

    out = ROOT / 'results' / ('type4_plugin_mock_head_%s.png' % head)
    fig.savefig(out, dpi=160, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print('saved:', out)
    if head == 'semiellipsoidal':
        (ROOT / 'results' / 'type4_plugin_dialog_mock.png').write_bytes(out.read_bytes())
        print('updated: type4_plugin_dialog_mock.png')


for h in HEADS:
    render(h)
