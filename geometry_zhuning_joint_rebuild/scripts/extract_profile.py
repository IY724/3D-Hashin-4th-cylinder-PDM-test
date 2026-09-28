# -*- coding: utf-8 -*-
"""从 zhuning_joint_rebuilt.cae 重建 θ=0 子午面轮廓（线/圆弧精确，样条按已知椭圆族重生），
链合成闭环，输出曲线定义与采样折线，供插件示意图与参数化内核使用。
运行于 geometry_zhuning_joint_rebuild 目录：
    E:\\ABAQUS2025\\Commands\\abaqus.BAT cae noGUI=scripts/extract_profile.py
输出：results/profile_extract.json（含体积校验）
"""
import json
import math
import os

from builtins import sum, min as _bmin
from abaqus import mdb

CAE = 'zhuning_joint_rebuilt.cae'
MODEL = 'Zhuning_Joint_Rebuilt'
TOL = 1e-6

A_E, B_E = 140.0, 70.0          # 外椭圆（paper-fit v3）
Y_C = 641.0
T_LINER = 5.0

VOL_REF = {'PA6_Liner': 193397.07, 'Al6061_Bosses': 60577.93}
SECTOR_RAD = math.radians(10.0)


def ellipse_point(u, off):
    x, y = A_E * math.cos(u), Y_C + B_E * math.sin(u)
    nx, ny = math.cos(u) / A_E, math.sin(u) / B_E
    nn = math.hypot(nx, ny)
    return (x + off * nx / nn, y + off * ny / nn)


def ellipse_param(r, y, off):
    """(r,y) 在偏置椭圆上求参数：粗扫 + 黄金分割细化；不在线上返回 None。"""
    def d(u):
        p = ellipse_point(u, off)
        return math.hypot(p[0] - r, p[1] - y)

    n = 4000
    lo, hi = -math.pi / 2, math.pi / 2
    best_u, best_d = lo, 1e18
    for k in range(n + 1):
        u = lo + (hi - lo) * k / float(n)
        du = d(u)
        if du < best_d:
            best_u, best_d = u, du
    a, b = best_u - (hi - lo) / n, best_u + (hi - lo) / n
    for _ in range(60):
        m1 = a + (b - a) * 0.382
        m2 = a + (b - a) * 0.618
        if d(m1) < d(m2):
            b = m2
        else:
            a = m1
    u = (a + b) / 2
    return u if d(u) < 1e-4 else None


def ellipse_pts(u1, u2, off, n=60):
    return [ellipse_point(u1 + (u2 - u1) * k / float(n), off) for k in range(n + 1)]


def dist(p, q):
    return math.hypot(p[0] - q[0], p[1] - q[1])


def circumcenter(A, B, P):
    ax, ay = A
    bx, by = B
    px, py = P
    d = 2.0 * (ax * (by - py) + bx * (py - ay) + px * (ay - by))
    if abs(d) < 1e-9:
        return None
    ux = ((ax * ax + ay * ay) * (by - py) + (bx * bx + by * by) * (py - ay) +
          (px * px + py * py) * (ay - by)) / d
    uy = ((ax * ax + ay * ay) * (px - bx) + (bx * bx + by * by) * (ax - px) +
          (px * px + py * py) * (bx - ax)) / d
    return (ux, uy)


def on_segment(A, B, P, tol=1e-6):
    return abs(dist(A, P) + dist(P, B) - dist(A, B)) < tol


def norm_angle(t):
    return (t + 2 * math.pi) % (2 * math.pi)


openMdb(pathName=CAE)
model = mdb.models[MODEL]
result = {'meta': {'cae': CAE, 'model': MODEL,
                   'note': 'theta=0 meridian; arcs exact via 3-point circle, '
                           'splines regenerated from known offset-ellipse families'}}

for part_name, key in (('PA6_Liner', 'liner'), ('Al6061_Bosses', 'boss')):
    part = model.parts[part_name]
    segments = []
    for edge in part.edges:
        vi = edge.getVertices()
        if len(vi) != 2:
            continue
        A3 = part.vertices[vi[0]].pointOn[0]
        B3 = part.vertices[vi[1]].pointOn[0]
        P3 = edge.pointOn[0]
        if abs(A3[2]) > 1e-7 or abs(B3[2]) > 1e-7 or abs(P3[2]) > 1e-7:
            continue          # 周向/斜面连接边，非子午线
        A = (A3[0], A3[1])
        B = (B3[0], B3[1])
        Pxy = (P3[0], P3[1])
        seg = None
        if on_segment(A, B, Pxy):
            seg = {'kind': 'line', 'A': A, 'B': B}
        else:
            cc = circumcenter(A, B, Pxy)
            if cc is not None:
                R = dist(cc, A)
                if 1e-3 < R < 1e5 and abs(dist(cc, B) - R) < 1e-4 and \
                        abs(dist(cc, Pxy) - R) < 1e-4:
                    a1 = math.atan2(A[1] - cc[1], A[0] - cc[0])
                    am = math.atan2(Pxy[1] - cc[1], Pxy[0] - cc[0])
                    a2 = math.atan2(B[1] - cc[1], B[0] - cc[0])
                    dm = norm_angle(am - a1)
                    db = norm_angle(a2 - a1)
                    if dm < db:
                        u1, u2 = a1, a1 + db
                    else:
                        u1, u2 = a1, a1 - (2 * math.pi - db)
                    curv = None
                    try:
                        cA = edge.getCurvature(point=A3)
                        cB = edge.getCurvature(point=B3)
                        kA = cA['curvature'] if isinstance(cA, dict) else cA
                        kB = cB['curvature'] if isinstance(cB, dict) else cB
                        curv = (kA, kB)
                    except Exception:
                        pass
                    seg = {'kind': 'arc', 'A': A, 'B': B, 'center': list(cc),
                           'radius': R, 'u1': u1, 'u2': u2, 'end_curvature': curv}
        if seg is None:
            fam = None
            for off in (0.0, -T_LINER, T_LINER):
                uA = ellipse_param(A[0], A[1], off)
                uB = ellipse_param(B[0], B[1], off)
                if uA is None or uB is None:
                    continue
                pm = ellipse_pts(uA, uB, off, 40)[20]
                if dist(pm, Pxy) < 0.05:
                    fam = (off, uA, uB)
                    break
            assert fam is not None, u'样条族匹配失败: %s A=%s B=%s' % (part_name, A, B)
            seg = {'kind': 'spline', 'A': A, 'B': B,
                   'ellipse_offset': fam[0], 'u1': fam[1], 'u2': fam[2]}
        segments.append(seg)

    # 端点链合成闭环
    pool = list(segments)
    loops = []
    while pool:
        s = pool.pop(0)
        chain = [s]
        start, end = s['A'], s['B']
        while dist(start, end) > TOL:
            hit = None
            for j, t in enumerate(pool):
                if dist(t['A'], end) < TOL:
                    hit = (j, t, t['B'])
                    break
                if dist(t['B'], end) < TOL:
                    flip = dict(t)
                    flip['A'], flip['B'] = t['B'], t['A']
                    if flip['kind'] == 'arc':
                        flip['u1'], flip['u2'] = t['u2'], t['u1']
                    hit = (j, flip, flip['B'])
                    break
            assert hit is not None, u'链中断于 %s (%s)' % (end, part_name)
            j, t, end = hit
            pool.pop(j)
            chain.append(t)
        loops.append(chain)

    def sample(seg):
        if seg['kind'] == 'line':
            return [seg['A'], seg['B']]
        if seg['kind'] == 'arc':
            u1, u2 = seg['u1'], seg['u2']
            n = max(8, int(abs(u2 - u1) * seg['radius'] / 0.4))
            return [(seg['center'][0] + seg['radius'] * math.cos(u1 + (u2 - u1) * k / float(n)),
                     seg['center'][1] + seg['radius'] * math.sin(u1 + (u2 - u1) * k / float(n)))
                    for k in range(n + 1)]
        return ellipse_pts(seg['u1'], seg['u2'], seg['ellipse_offset'])

    poly_loops = []
    for chain in loops:
        pts = []
        for seg in chain:
            sp = sample(seg)
            pts += sp if not pts else sp[1:]
        dedup = [pts[0]]
        for q in pts[1:]:
            if dist(q, dedup[-1]) > 1e-9:
                dedup.append(q)
        poly_loops.append(dedup)

    moment = 0.0
    for pts in poly_loops:
        m = 0.0
        for k in range(len(pts)):
            cr = pts[k][0] * pts[(k + 1) % len(pts)][1] - pts[(k + 1) % len(pts)][0] * pts[k][1]
            m += (pts[k][0] + pts[(k + 1) % len(pts)][0]) * cr / 6.0
        moment += abs(m)          # 两阀座环路绕向相反，逐环取绝对值
    volume = moment * SECTOR_RAD
    ref = VOL_REF[part_name]
    vol_err = abs(volume - ref) / ref
    result[key] = {'segments': segments, 'loops': poly_loops,
                   'n_loops': len(loops), 'volume_mm3': volume,
                   'volume_ref_mm3': ref, 'volume_rel_err': vol_err}
    print('%s: %d segs, %d loops, volume %.2f vs ref %.2f (err %.2e)' %
          (part_name, len(segments), len(loops), volume, ref, vol_err))

ok = all(d['volume_rel_err'] < 1e-3 for k, d in result.items() if k in VOL_REF)
result['meta']['volume_check_passed'] = bool(ok)
os.makedirs('results', exist_ok=True)
with open(os.path.join('results', 'profile_extract.json'), 'w', encoding='utf-8') as f:
    json.dump(result, f, ensure_ascii=False, indent=1)
print('PROFILE_EXTRACT_DONE check=%s' % ok)
