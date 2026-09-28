# -*- coding: utf-8 -*-
"""Type-IV 氢气瓶建模插件 · 对话框设计样图 v3（原创两栏版式）：
左栏=封头类型选择（半球/半椭球/等张力）+ 分组参数填空（含角度输入）；
右栏=符号示意图（随输入自动重画）。无 Model/Part 输入：插件在当前模型内
自动创建并固定命名 Liner / Boss 两个 Part。不需要 Abaqus。

运行：python scripts/draw_plugin_dialog_mock.py
输出：results/type4_plugin_dialog_mock.png
"""
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch, Circle

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

ROOT = Path(__file__).resolve().parent.parent
DATA = json.loads((ROOT / 'profiles.json').read_text(encoding='utf-8'))

C = '#1a3fbf'
CW = '#263442'
UI = '#263442'
BOX = '#8b99a6'
GRAY = '#8a949e'

fig, ax = plt.subplots(figsize=(11.6, 7.2))
ax.set_xlim(0, 100)
ax.set_ylim(0, 61)
ax.set_aspect('equal')
ax.axis('off')

ax.add_patch(Rectangle((0.5, 0.5), 99, 60, fc='white', ec=BOX, lw=1.2))
ax.add_patch(Rectangle((0.5, 57.6), 99, 2.9, fc='#dfe5ea', ec=BOX, lw=1.2))
ax.text(2.2, 59.05, 'Type-IV Hydrogen Vessel  —  Liner & Boss Generator', fontsize=11,
        weight='bold', color=UI, va='center')

# ================= 左栏：参数输入 =================
ax.text(3.0, 51.9, 'Inputs', fontsize=9.5, weight='bold', color=UI)
ax.add_patch(Rectangle((2, 19.5), 38, 30.6, fc='#f4f6f8', ec=BOX, lw=1.0))

def group(y, name):
    ax.text(3.2, y, name, fontsize=8.5, weight='bold', color=UI, va='center')
    ax.plot([3.2, 38.8], [y - 0.78] * 2, color='#c3ccd4', lw=0.8)

def mfield(x, y, sym, val, lw=3.0, w=8.2, off=False):
    col = GRAY if off else UI
    ax.text(x, y, '$%s$' % sym, fontsize=8.5, color=col, va='center')
    ax.add_patch(Rectangle((x + lw, y - 0.95), w, 2.0,
                           fc='#f0f0f0' if off else 'white', ec=BOX, lw=0.9))
    ax.text(x + lw + w - 0.45, y, val, fontsize=8, color=GRAY if off else '#555f6a',
            va='center', ha='right')

# ---- 封头类型（三种） ----
group(49.2, 'Head type')
heads = [('Hemispherical', lambda u: 1.35 * np.sqrt(np.maximum(0, 1 - u ** 2)), False),
         ('Semi-ellipsoidal', lambda u: 0.95 * np.sqrt(np.maximum(0, 1 - u ** 2)), True),
         ('Isotensoid', lambda u: 1.15 * np.maximum(0, 1 - u ** 4) ** 0.38, False)]
for k, (name, f, sel) in enumerate(heads):
    yc = 47.4 - k * 1.62
    ax.add_patch(Circle((3.7, yc + 0.5), 0.42, fc=C if sel else 'white', ec=BOX, lw=0.9))
    u = np.linspace(-1, 1, 30)
    ax.plot(4.9 + 1.5 * u, yc + 0.35 + f(u), color=C if sel else GRAY, lw=1.2)
    ax.text(8.1, yc + 0.55, name, fontsize=8, color=C if sel else '#5a646e',
            va='center', weight='bold' if sel else 'normal')
mfield(3.4, 42.3, 'k', '2.0')
mfield(20.0, 42.3, '\\alpha_0', '55', off=True)
ax.text(3.4, 40.6, 'fields switch with head type  (Isotensoid uses $r_0$ below)',
        fontsize=6.8, color='#5a646e')

# ---- 分组参数 ----
group(38.9, 'Barrel  (cylindrical section)')
for i, (s, v) in enumerate([('R', '135'), ('l_1', '1282'), ('t_1', '5')]):
    mfield(3.4 + i * 11.4, 36.6, s, v)

group(34.6, 'Shoulder / dome transition')
for i, (s, v) in enumerate([('l_2', '119'), ('t_2', '5'), ('r_1', '10')]):
    mfield(3.4 + i * 11.4, 32.3, s, v)

group(30.3, 'Polar joint  (liner-boss bond)')
for i, (s, v) in enumerate([('r_f', '46'), ('l_c', '40')]):
    mfield(3.4 + i * 16.6, 28.0, s, v, lw=3.4, w=9.6)
for i, (s, v) in enumerate([('t_3', '5'), ('r_2', '10')]):
    mfield(3.4 + i * 16.6, 25.9, s, v, lw=3.4, w=9.6)

group(23.9, 'Polar boss  (metal insert)')
for i, (s, v) in enumerate([('r_0', '16'), ('r_b', '36'), ('t_4', '20')]):
    mfield(3.4 + i * 11.4, 21.6, s, v)

# ================= 右栏：示意图 =================
ax.text(43.5, 51.9, 'Schematic  (redrawn from inputs)', fontsize=9.5, weight='bold', color=UI)
ax.add_patch(Rectangle((42, 19.5), 56, 30.6, fc='#ebebeb', ec=BOX, lw=1.0))

X = lambda x: 43.5 + 0.163 * x     # x = 距极孔端面轴向距离 mm
Y = lambda r: 22.5 + 0.150 * r     # r = 半径 mm
end = 760.0
CUT = 300.0


def loop_pts(loop):
    return [tuple(p) for c in loop for p in c['points'][:-1]] + [tuple(loop[-1]['points'][-1])]


def cut_and_transform(loop, ymin):
    pts = loop_pts(loop)
    n = len(pts)
    crossings = []
    for i in range(n):
        y1, y2 = pts[i][1], pts[(i + 1) % n][1]
        if (y1 - ymin) * (y2 - ymin) < 0:
            t = (ymin - y1) / (y2 - y1)
            r = pts[i][0] + t * (pts[(i + 1) % n][0] - pts[i][0])
            crossings.append((i, (r, ymin)))
    assert len(crossings) == 2
    (i1, c1), (i2, c2) = crossings
    chain = [c1]
    k = (i1 + 1) % n
    while True:
        chain.append(pts[k])
        if k == i2:
            break
        k = (k + 1) % n
    chain.append(c2)
    if max(p[1] for p in chain) < ymin + 50:
        chain = [c2]
        k = (i2 + 1) % n
        while True:
            chain.append(pts[k])
            if k == i1:
                break
            k = (k + 1) % n
        chain.append(c1)
    return np.array([[X(end - y), Y(r)] for r, y in chain])


liner_P = cut_and_transform(DATA['liner'], 430.0)
boss_P = np.array([[X(end - y), Y(r)] for r, y in loop_pts(DATA['boss_positive'])])
for P, fc in ((liner_P, '#cfe0f2'), (boss_P, '#f6dfc2')):
    ax.add_patch(plt.Polygon(P, closed=True, fc=fc, ec=CW, lw=1.5, zorder=2))

ax.plot([42.6, 97.4], [Y(0)] * 2, color='#e08bb0', lw=1.0, ls=(0, (7, 3, 1, 3)), zorder=1)
zx = X(CUT + 10)
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

dim_h(X(0), X(119), 47.4, '$l_2$', ty=48.0)
dim_h(X(119), X(CUT), 49.1, '$l_1$', ty=49.7)
dim_v(Y(0), Y(135), 95.0, '$R$')
dim_h(X(0), X(CUT), 21.2, '$L$', ty=20.4)

leader(X(30), Y(16), 44.2, 23.8, '$r_0$')
leader(X(30), Y(26), 45.2, 29.9, '$t_4$')
leader(X(30), Y(36), 46.5, 33.2, '$r_b$')
leader(X(64), Y(44), 50.0, 42.7, '$r_1$')
leader(X(90), Y(124), 54.5, 44.9, '$t_2$')
leader(X(112), Y(61), 73.0, 33.6, '$l_c$')
leader(X(125), Y(49), 75.0, 30.9, '$t_3$')
leader(X(135), Y(46), 76.0, 28.3, '$r_f$')
leader(X(114), Y(21), 77.0, 25.4, '$r_2$')
leader(X(250), Y(137.5), 90.0, 45.9, '$t_1$')

# ================= 底部 =================
ax.text(2.5, 16.6, 'Sector angle  $\\theta$:', fontsize=8.5, color=UI, va='center')
ax.add_patch(Rectangle((12.8, 15.65), 7.0, 2.0, fc='white', ec=BOX, lw=0.9))
ax.text(19.35, 16.65, '10', fontsize=8, color='#555f6a', va='center', ha='right')
ax.text(20.8, 16.65, '°   (360 = full model)', fontsize=8, color='#5a646e', va='center')
ax.text(30.5, 16.65, 'Reference standard:  GB/T 42612-2023', fontsize=8.5, color=UI, va='center')
ax.text(2.5, 13.6, "Parts are created in the current model and named automatically:  'Liner' ,  'Boss'.\n"
                   "All dimensions in mm; values shown are an editable preset.",
        fontsize=7.5, color='#5a646e', va='center')

for x0, w, name in [(68, 9, 'OK'), (78.5, 11, 'Apply'), (91, 8, 'Cancel')]:
    ax.add_patch(FancyBboxPatch((x0, 1.4), w, 2.9, boxstyle='round,pad=0.25',
                                fc='#e4e9ee', ec=BOX, lw=1.0))
    ax.text(x0 + w / 2, 2.85, name, fontsize=9, color=UI, ha='center', va='center')

out = ROOT / 'results' / 'type4_plugin_dialog_mock.png'
fig.savefig(out, dpi=160, bbox_inches='tight', facecolor='white')
print('saved:', out)
