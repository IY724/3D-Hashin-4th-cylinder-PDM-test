# -*- coding: utf-8 -*-
"""Ⅳ型储氢气瓶（内胆+瓶阀）参数标注图：按 profiles.json 真实轮廓绘制上半剖面
（布局参照常规Ⅳ型瓶内胆图纸：极孔/折边段/肩部/筒身逐项标注），
标出插件需要用户填写的尺寸参数（单位 mm）。不需要 Abaqus。

运行：python scripts/draw_parameter_drawing.py
输出：results/type4_parameter_drawing.png
"""
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parent.parent
DATA = json.loads((ROOT / 'profiles.json').read_text(encoding='utf-8'))
P = DATA['parameters']
d, j = P['paper_dimensions'], P['fitted_joint_dimensions']

a = DATA['derived']['ellipse_a']
b = DATA['derived']['ellipse_b']
h = d['cylinder_length'] / 2.0
end = d['overall_length'] / 2.0
t = d['liner_nominal_thickness']
neck = d['pole_outer_diameter'] / 2.0
bore = d['bore_diameter'] / 2.0
rf = j['neck_outer_fillet_radius']
rb = j['bore_inner_fillet_radius']
ps = DATA['derived']['split_point']
liner_half = ps[1]
collar_tip = h + j['collar_tip_from_tangent']          # 621
flare_top = h + j['collar_return_center_from_tangent'] + j['collar_return_radius']  # 661
rc = j['collar_interface_radius']

C = '#1a3fbf'
CW = '#263442'


def poly(loop):
    return np.array([p for c in loop for p in c['points'][:-1]] + [loop[-1]['points'][-1]])


liner_pts = poly(DATA['liner'])
boss_pts = poly(DATA['boss_positive'])
boss_neg_pts = poly(DATA['boss_negative'])


def ellipse_pt(u, offset=0.0):
    x, y = a * math.cos(u), h + b * math.sin(u)
    nx, ny = math.cos(u) / a, math.sin(u) / b
    n = math.hypot(nx, ny)
    return x + offset * nx / n, y + offset * ny / n


plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
fig, ax = plt.subplots(figsize=(7.6, 12.8))

ax.fill(liner_pts[:, 0], liner_pts[:, 1], color='#cfe0f2', alpha=.85, zorder=1)
ax.fill(boss_pts[:, 0], boss_pts[:, 1], color='#f6dfc2', alpha=.9, zorder=1)
ax.fill(boss_neg_pts[:, 0], boss_neg_pts[:, 1], color='#f6dfc2', alpha=.9, zorder=1)
ax.plot(liner_pts[:, 0], liner_pts[:, 1], color=CW, lw=2.0, zorder=3)
ax.plot(boss_pts[:, 0], boss_pts[:, 1], color=CW, lw=1.3, zorder=3)
ax.plot(boss_neg_pts[:, 0], boss_neg_pts[:, 1], color=CW, lw=1.3, zorder=3)

ax.axvline(0, color='#e08bb0', lw=1.0, ls=(0, (7, 3, 1, 3)), zorder=2)
ax.plot([-14, 150], [0, 0], color='#e08bb0', lw=1.0, ls=(0, (7, 3, 1, 3)), zorder=2)
ax.text(58, 16, '中面（对称）', fontsize=8.5, color='#c2507e', ha='left', va='bottom')

DIM = dict(arrowstyle='<|-|>', color=C, lw=0.95, mutation_scale=9, shrinkA=0, shrinkB=0)


def dim_h(x1, x2, y, ext=None, text='', tx=None, ty=None, ta='center'):
    for xe, yf in (ext or []):
        ax.plot([xe, xe], [min(yf, y) - 2, max(yf, y) + 2], lw=0.55, color=C, zorder=2)
    ax.annotate('', xy=(x1, y), xytext=(x2, y), arrowprops=DIM, zorder=4)
    ax.text(tx if tx is not None else (x1 + x2) / 2, ty if ty is not None else y, text,
            ha=ta, va='center', color=C, fontsize=10.5, zorder=5)


def dim_v(y1, y2, x, ext=None, text='', tx=None):
    for ye, xf in (ext or []):
        ax.plot([min(xf, x) - 2, max(xf, x) + 2], [ye, ye], lw=0.55, color=C, zorder=2)
    ax.annotate('', xy=(x, y1), xytext=(x, y2), arrowprops=DIM, zorder=4)
    ax.text(tx if tx is not None else x + 5, (y1 + y2) / 2, text, ha='left', va='center',
            color=C, fontsize=10.5, rotation=90, zorder=5)


def leader(xy, xytext, text, ha='left'):
    ax.annotate(text, xy=xy, xytext=xytext, fontsize=10.5, color=C, ha=ha, va='center',
                arrowprops=dict(arrowstyle='->', color=C, lw=0.8, shrinkA=2, shrinkB=1),
                zorder=5)


# ---- 极孔/瓶口（顶部） ----
dim_h(0, neck, 776, ext=[(neck, end)], text='瓶口外径 Ø%.0f' % (2 * neck), tx=42, ta='left')
dim_h(0, bore, 806, ext=[(bore, end)], text='极孔内径 Ø%.0f' % (2 * bore), tx=20, ta='left')

pf = DATA['derived']['ellipse_fillet_point']
fc = DATA['derived']['fillet_center']
leader(((pf[0] + fc[0]) / 2, (pf[1] + fc[1]) / 2), (-12, 742), '极孔外圆角 R%.0f' % rf)
leader((bore + rb * (1 - 0.7071), 641 + rb * (1 - 0.7071)), (-12, 698), '极孔内圆角 R%.0f' % rb)

dim_h(0, rc, 585, ext=[(rc, collar_tip)], text='极孔外径 Ø%.0f *' % (2 * rc), tx=rc / 2, ty=593)

# ---- 肩部（椭球封头）法向壁厚 ----
u_sh = 0.45
ox, oy = ellipse_pt(u_sh)
ix, iy = ellipse_pt(u_sh, -t)
ax.annotate('', xy=(ix, iy), xytext=(ox, oy), arrowprops=DIM, zorder=4)
leader(((ix + ox) / 2, (iy + oy) / 2), (166, 710), '肩部厚度 %.0f' % t)

# ---- 折边段（极部回折段） ----
ax.annotate('', xy=(46.6, 624.5), xytext=(50.6, 624.5), arrowprops=DIM, zorder=4)
leader((48.6, 624.5), (58, 556), '折边段厚度 %.0f *' % t)
dim_v(collar_tip, flare_top, 104, ext=[(rc + t, collar_tip), (93, flare_top)],
      text='折边段长度 %.0f *' % (flare_top - collar_tip), tx=109)

# ---- 筒身 ----
ax.annotate('', xy=(135, -180), xytext=(140, -180), arrowprops=DIM, zorder=4)
leader((137.5, -180), (154, -216), '筒身壁厚 %.0f' % t)

# ---- 底部 ----
dim_h(0, a, -28, ext=[(a, -44)], text='内胆外径 Ø%.0f' % (2 * a), tx=a / 2, ty=-41)

# ---- 图例 ----
ax.legend(handles=[Patch(fc='#cfe0f2', ec=CW, label='PA6 内胆'),
                   Patch(fc='#f6dfc2', ec=CW, label='Al6061 阀座')],
          loc='lower left', bbox_to_anchor=(0.02, 0.115), fontsize=9.5, framealpha=.92)

ax.set_title('Ⅳ型储氢气瓶参数标注图（内胆 + 瓶阀）—— 插件输入参数', fontsize=13.5,
             color=CW, pad=10)
ax.set_xlabel('半径 r / mm', fontsize=10)
ax.set_xlim(-14, 248)
ax.set_ylim(-46, 836)
ax.set_aspect('equal')
ax.set_yticks([])
ax.grid(False)
for s in ax.spines.values():
    s.set_visible(False)
ax.tick_params(labelsize=9)

note = ('封头类型（插件下拉项）：椭球形 a/b=%.1f（本图）｜半球形｜碟形｜等张力 isotensoid；扇区角（填空）：10°（可填 360° 整周）\n'
        '全长（对称）：筒身段 %.0f｜内胆 %.0f *｜整瓶 %.0f；瓶口外露 ≈ %.0f *；单位 mm；带 * 为拟合值（论文未标注，可修改）' % (
            d['ellipse_axis_ratio'], 2 * h, 2 * liner_half, 2 * end, end - liner_half))
fig.text(0.5, 0.004, note, fontsize=9.2, color=CW, ha='center', va='bottom')

fig.savefig(ROOT / 'results' / 'type4_parameter_drawing.png', dpi=170, bbox_inches='tight')
print('saved:', ROOT / 'results' / 'type4_parameter_drawing.png')
