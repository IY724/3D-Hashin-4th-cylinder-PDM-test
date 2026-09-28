# -*- coding: utf-8 -*-
"""生成可复现截面、尺寸来源和几何核对图；不需要 Abaqus。"""
from pathlib import Path
import copy
import json
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent  # 脚本位于 scripts\，数据根目录为上级
P = json.loads((ROOT / 'parameters.json').read_text(encoding='utf-8'))
d, j = P['paper_dimensions'], P['fitted_joint_dimensions']
a = d['liner_outer_diameter'] / 2
b = a / d['ellipse_axis_ratio']
h = d['cylinder_length'] / 2
end = d['overall_length'] / 2
t = d['liner_nominal_thickness']
neck = d['pole_outer_diameter'] / 2
bore = d['bore_diameter'] / 2
rs = j['boss_flange_split_radius']
rc = j['collar_interface_radius']
rcurve = j['collar_return_radius']
cx = rc + rcurve
cy = h + j['collar_return_center_from_tangent']
lip = h + j['collar_tip_from_tangent']
fy = cy + rcurve
rf = j['neck_outer_fillet_radius']
rb = j['bore_inner_fillet_radius']


def root(function, lo=0.0, hi=math.pi / 2):
    flo = function(lo)
    assert flo * function(hi) <= 0, '根未被区间包围'
    for unused in range(80):
        mid = (lo + hi) / 2
        if flo * function(mid) <= 0:
            hi = mid
        else:
            lo, flo = mid, function(mid)
    return (lo + hi) / 2


def ellipse(u, offset=0):
    x, y = a * math.cos(u), h + b * math.sin(u)
    nx, ny = math.cos(u) / a, math.sin(u) / b
    norm = math.hypot(nx, ny)
    return [x + offset * nx / norm, y + offset * ny / norm]


def line(p0, p1):
    assert math.dist(p0, p1) > 1e-8
    return {'kind': 'line', 'points': [list(p0), list(p1)]}


def arc(center, radius, start, stop):
    values = np.linspace(math.radians(start), math.radians(stop), 81)
    points = [[center[0] + radius * math.cos(u), center[1] + radius * math.sin(u)] for u in values]
    return {'kind': 'arc', 'center': list(center), 'points': points,
            'ccw': stop > start, 'radius': radius}


def spline(u0, u1, offset=0):
    return {'kind': 'spline', 'points': [ellipse(u, offset) for u in np.linspace(u0, u1, 121)]}


def reverse(curves):
    result = copy.deepcopy(list(reversed(curves)))
    for c in result:
        c['points'].reverse()
        if c['kind'] == 'arc':
            c['ccw'] = not c['ccw']
    return result


def reflect(curves):
    result = copy.deepcopy(curves)
    for c in result:
        c['points'] = [[x, -y] for x, y in c['points']]
        if c['kind'] == 'arc':
            c['center'][1] *= -1
            c['ccw'] = not c['ccw']
    return result


us = math.acos(rs / a)
uf = root(lambda u: ellipse(u, rf)[0] - (neck + rf))
fillet_center = ellipse(uf, rf)
pf = ellipse(uf)
fillet_angle = math.degrees(math.atan2(pf[1] - fillet_center[1], pf[0] - fillet_center[0]))
blend_r = j['liner_inner_blend_radius']
ui = root(lambda u: ellipse(u, -t - blend_r)[0] - (rs + t + blend_r))
blend_center = ellipse(ui, -t - blend_r)
ps = ellipse(us)
pi = ellipse(ui, -t)
blend_angle = math.degrees(math.atan2(pi[1] - blend_center[1], pi[0] - blend_center[0]))
assert t < rcurve and rs > cx and fy < blend_center[1] < pi[1] < ps[1]
assert fillet_center[1] < end and neck > bore

# PA6 与铝阀座严格复用的共同轮廓。
interface = [line(ps, (rs, fy)), line((rs, fy), (cx, fy)),
             arc((cx, cy), rcurve, 90, 180), line((rc, cy), (rc, lip))]
positive_liner = [spline(0, us)] + interface + [
    line((rc, lip), (rc + t, lip)), line((rc + t, lip), (rc + t, cy)),
    arc((cx, cy), rcurve - t, 180, 90),
    line((cx, fy - t), (rs, fy - t)), arc((rs, fy), t, -90, 0),
    line((rs + t, fy), (rs + t, blend_center[1])),
    arc(blend_center, blend_r, 180, blend_angle), spline(ui, 0, -t)]
liner = [line((a, -h), (a, h))] + positive_liner + [line((a - t, h), (a - t, -h))] + reverse(reflect(positive_liner))

outer_boss = [spline(us, uf), arc(fillet_center, rf, fillet_angle, -180)]
boss = outer_boss + [line((neck, fillet_center[1]), (neck, end)),
    line((neck, end), (bore, end)), line((bore, end), (bore, h + rb)),
    arc((bore + rb, h + rb), rb, 180, 270),
    line((bore + rb, h), (neck, h)), line((neck, h), (neck, lip)),
    line((neck, lip), (rc, lip))] + reverse(interface)
outer_positive = [spline(0, us)] + outer_boss
wcm = reverse(reflect(outer_positive)) + [line((a, -h), (a, h))] + outer_positive


def sample(curves, closed=True):
    for c0, c1 in zip(curves, curves[1:] + (curves[:1] if closed else [])):
        assert math.dist(c0['points'][-1], c1['points'][0]) < 1e-7, '截面轮廓不连续'
    return np.array([p for c in curves for p in c['points'][:-1]] + [curves[-1]['points'][-1]])


def measure(curves):
    xy = sample(curves)
    x, y = xy[:, 0], xy[:, 1]
    cr = x[:-1] * y[1:] - x[1:] * y[:-1]
    area = abs(float(np.sum(cr) / 2))
    moment = abs(float(np.sum((x[:-1] + x[1:]) * cr) / 6))
    # 排除非相邻线段真正相交；相邻端点属于合法连接。
    for i in range(len(xy) - 2):
        p, q = xy[i], xy[i + 1]
        u, v = xy[i + 2:-1], xy[i + 3:]
        if i == 0:
            u, v = u[:-1], v[:-1]
        cross = lambda aa, bb: aa[..., 0] * bb[..., 1] - aa[..., 1] * bb[..., 0]
        s1, s2 = cross(q - p, u - p), cross(q - p, v - p)
        s3, s4 = cross(v - u, p - u), cross(v - u, q - u)
        assert not np.any((s1 * s2 < -1e-12) & (s3 * s4 < -1e-12)), '截面自相交'
    return {'sampled_area_mm2': area,
            'sampled_sector_volume_mm3': moment * math.radians(d['sector_angle_deg']),
            'closed': True, 'proper_self_intersections': 0}


profiles = {'parameters': P, 'derived': {'ellipse_a': a, 'ellipse_b': b,
    'tangent_y': h, 'fillet_center': fillet_center, 'ellipse_fillet_point': pf,
    'split_point': ps, 'inner_trim_point': pi, 'inner_blend_center': blend_center, 'collar_center': [cx, cy]},
    'liner': liner, 'boss_positive': boss, 'boss_negative': reflect(boss),
    'wcm_profile': wcm, 'joint_interface_positive': interface,
    'checks': {'liner': measure(liner), 'boss_positive': measure(boss)}}
profiles['checks']['nominal_thickness_check_mm'] = [math.dist(ellipse(u), ellipse(u, -t)) for u in np.linspace(0, ui, 17)]
(ROOT / 'profiles.json').write_text(json.dumps(profiles, ensure_ascii=False, indent=2), encoding='utf-8')
with (ROOT / 'wcm_meridian.csv').open('w', encoding='utf-8') as f:
    f.write('radius_mm,axial_y_mm\n')
    for r, y in sample(wcm, closed=False):
        f.write('{:.9f},{:.9f}\n'.format(r, y))

# 几何核对图，不属于计算结果或论文应力图。
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
fig, axes = plt.subplots(2, 1, figsize=(15, 9), gridspec_kw={'height_ratios': [1, 2]})
colors = ['#3287b0', '#d89847', '#d89847']
for ax in axes:
    for curves, color, label in zip([liner, boss, reflect(boss)], colors, ['PA6 内衬', '6061-T6 铝阀座', None]):
        pts = sample(curves)
        ax.fill(pts[:, 1], pts[:, 0], color=color, alpha=.75, label=label)
        ax.plot(pts[:, 1], pts[:, 0], color='#263442', lw=.8)
    ax.set_aspect('equal')
    ax.set_xlabel('轴向 Y / mm')
    ax.set_ylabel('半径 r / mm')
    ax.grid(alpha=.2)
axes[0].set_xlim(-780, 780)
axes[0].set_ylim(0, 155)
axes[0].set_title('朱宁基体截面：全长 1520，筒身 1282，外径 280，旋转 10°')
axes[0].legend(loc='lower center', ncol=2)
axes[1].set_xlim(595, 775)
axes[1].set_ylim(0, 148)
axes[1].axvline(h, color='#888', ls='--', lw=.8)
axes[1].set_title('端部放大：材料界面共用；橙色与蓝色不重叠。接头未标注尺寸为拟合值。')
fig.tight_layout()
fig.savefig(ROOT / 'section_check.png', dpi=180)
plt.close(fig)
print(json.dumps({'derived': profiles['derived'], 'checks': profiles['checks']}, ensure_ascii=False, indent=2))

# 网页返回的国标扫描本已在本地缓存，仅抽取核对页；不可据此声称规范认证。
standard_cache = Path(r'C:\Users\Administrator\.qoder\cache\projects\3D-Hashin-4th-cylinder-PDM-8b93d1eb\agent-tools\e69740e5\8f649272.txt')
if standard_cache.exists():
    import fitz
    with fitz.open(str(standard_cache), filetype='pdf') as doc:
        for number in (1, 8, 9, 10, 11, 12, 13, 14):
            doc[number - 1].get_pixmap(matrix=fitz.Matrix(1.7, 1.7)).save(str(ROOT / 'evidence' / ('gbt42612_page_{}.png'.format(number))))
