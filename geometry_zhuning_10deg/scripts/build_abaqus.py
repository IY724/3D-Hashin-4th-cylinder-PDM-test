# -*- coding: utf-8 -*-
"""Abaqus/CAE 2025 原生建模：二维闭合截面旋转10°，保留两种材料。

运行前先用普通 Python 执行 prepare_geometry.py。
本脚本只创建新模型和带时间戳的输出目录，不提交分析作业。
"""
from abaqus import mdb, session
from abaqusConstants import *
import caeModules
from builtins import sum, min, max
import datetime
import json
import math
import os
import traceback

# Abaqus 的 noGUI 执行环境可能不提供 __file__，此时使用启动工作目录。
# 脚本位于 scripts\，数据根目录为其上级；若在工作目录直接运行请把 ROOT 改为当前目录。
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) if '__file__' in globals() else os.getcwd()
with open(os.path.join(ROOT, 'profiles.json'), 'r', encoding='utf-8') as f:
    DATA = json.load(f)
P = DATA['parameters']
D = P['paper_dimensions']
RUN = os.path.join(ROOT, 'runs', datetime.datetime.now().strftime('%Y%m%d_%H%M%S'))
os.makedirs(RUN)
os.chdir(RUN)
REPORT = {'status': 'BUILDING', 'run_directory': RUN, 'analysis_submitted': False,
          'parameters': P, 'warnings': []}


def write_report():
    with open(os.path.join(RUN, 'geometry_validation.json'), 'w', encoding='utf-8') as f:
        json.dump(REPORT, f, ensure_ascii=False, indent=2)


def draw(sketch, curves):
    for c in curves:
        points = tuple(tuple(p) for p in c['points'])
        if c['kind'] == 'line':
            sketch.Line(point1=points[0], point2=points[-1])
        elif c['kind'] == 'arc':
            sketch.ArcByCenterEnds(center=tuple(c['center']), point1=points[0],
                                  point2=points[-1], direction=COUNTERCLOCKWISE if c['ccw'] else CLOCKWISE)
        else:
            sketch.Spline(points=points, constrainPoints=False)


def sketch_for(model, name, loops):
    sketch = model.ConstrainedSketch(name=name, sheetSize=3200.0)
    sketch.ConstructionLine(point1=(0.0, -1000.0), point2=(0.0, 1000.0))
    for loop in loops:
        draw(sketch, loop)
    return sketch


def volume(part):
    # Cell.getSize 的快速体积估计不足以验证薄壁曲面；使用高精度质量属性积分。
    return part.getMassProperties(relativeAccuracy=HIGH)['volume']


def make_part(model, name, loops, material, set_name):
    sk = sketch_for(model, 'Sketch_' + name, loops)
    part = model.Part(name=name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
    part.BaseSolidRevolve(sketch=sk, angle=D['sector_angle_deg'], flipRevolveDirection=ON)
    region = part.Set(name=set_name, cells=part.cells)
    part.SectionAssignment(region=region, sectionName='Section_' + material)
    part.checkGeometry(detailed=ON, level=20)
    return part


def face_array(part, indices):
    result = part.faces[0:0]
    for i in indices:
        result += part.faces[i:i + 1]
    return result


def edge_array(part, indices):
    result = part.edges[0:0]
    for i in indices:
        result += part.edges[i:i + 1]
    return result


def on_wcm(xyz, tolerance=0.002):
    x, y, z = xyz
    r, y = math.hypot(x, z), abs(y)
    h = DATA['derived']['tangent_y']
    a, b = DATA['derived']['ellipse_a'], DATA['derived']['ellipse_b']
    fc = DATA['derived']['fillet_center']
    py = DATA['derived']['ellipse_fillet_point'][1]
    if y <= h + 1e-6:
        return abs(r - a) < tolerance
    if y <= py + 1e-6:
        return abs(r * r / a ** 2 + (y - h) ** 2 / b ** 2 - 1) < 2e-5
    if y <= fc[1] + 1e-6:
        return abs(math.hypot(r - fc[0], y - fc[1]) - P['fitted_joint_dimensions']['neck_outer_fillet_radius']) < tolerance
    return False


def main():
    name = P['model_name']
    if name in mdb.models:
        raise RuntimeError('目标模型已存在；请在新的 CAE 会话运行，避免覆盖已有工作。')
    model = mdb.Model(name=name)
    for key in ('PA6', 'Al6061_T6'):
        props = P['materials'][key]
        mat = model.Material(name=key)
        mat.Elastic(table=((props['young_modulus'], props['poisson_ratio']),))
        mat.Density(table=((props['density'],),))
        model.HomogeneousSolidSection(name='Section_' + key, material=key)

    liner = make_part(model, 'PA6_Liner', [DATA['liner']], 'PA6', 'MAT_PA6')
    bosses = make_part(model, 'Al6061_Bosses', [DATA['boss_positive'], DATA['boss_negative']], 'Al6061_T6', 'MAT_AL6061')
    assert len(liner.cells) == 1 and len(bosses.cells) == 2
    source_volumes = {'PA6': volume(liner), 'Al6061_T6': volume(bosses)}
    assembly = model.rootAssembly
    assembly.DatumCsysByDefault(CARTESIAN)
    i1 = assembly.Instance(name='PA6_Liner-1', part=liner, dependent=ON)
    i2 = assembly.Instance(name='Al6061_Bosses-1', part=bosses, dependent=ON)
    merged_instance = assembly.InstanceFromBooleanMerge(
        name='Vessel_Base', instances=(i1, i2), keepIntersections=ON,
        originalInstances=SUPPRESS, domain=GEOMETRY)
    part = model.parts[merged_instance.partName]
    part.checkGeometry(detailed=ON, level=20)
    assert len(part.cells) == 3, '合并后应为一个PA6实体区域和两个铝阀座实体区域'
    assert len(part.sets['MAT_PA6'].cells) == 1
    assert len(part.sets['MAT_AL6061'].cells) == 2
    # 显式覆盖合并自动映射，确保材料区不被合并成同一种材料。
    part.SectionAssignment(region=part.sets['MAT_PA6'], sectionName='Section_PA6')
    part.SectionAssignment(region=part.sets['MAT_AL6061'], sectionName='Section_Al6061_T6')

    coordinates = [v.pointOn[0] for v in part.vertices]
    angles = [math.degrees(math.atan2(z, x)) for x, y, z in coordinates]
    extent = [min(q[1] for q in coordinates), max(q[1] for q in coordinates)]
    radii = [math.hypot(q[0], q[2]) for q in coordinates]
    merged_volume = volume(part)
    loss = sum(source_volumes.values()) - merged_volume
    relative_loss = abs(loss) / sum(source_volumes.values())
    assert relative_loss < 1e-8, '合并发生体积损失，可能存在穿模'
    assert abs(extent[1] - extent[0] - D['overall_length']) < 1e-6
    assert abs(max(radii) * 2 - D['liner_outer_diameter']) < 1e-6
    assert abs(min(radii) * 2 - D['bore_diameter']) < 1e-6
    assert abs(min(angles) + D['sector_angle_deg']) < 1e-5 and abs(max(angles)) < 1e-5

    # 共享面建立实体邻接图；三个区域必须组成同一个连通体。
    graph = {c.index: set() for c in part.cells}
    internal, side0, side10, external, inside = [], [], [], [], []
    for face in part.faces:
        adjacent = tuple(face.getCells())
        if len(adjacent) == 2:
            internal.append(face.index)
            graph[adjacent[0]].add(adjacent[1])
            graph[adjacent[1]].add(adjacent[0])
            continue
        xyz = face.pointOn[0]
        x, y, z = xyz
        normal = face.getNormal(point=xyz)
        theta = math.atan2(z, x)
        circumferential = (-math.sin(theta), 0.0, math.cos(theta))
        plane_normal = abs(sum(a * b for a, b in zip(normal, circumferential)))
        if plane_normal > .99999:
            (side0 if abs(z) < 1e-6 else side10).append(face.index)
        elif on_wcm(xyz):
            external.append(face.index)
        else:
            r = math.hypot(x, z)
            is_neck = abs(r - D['pole_outer_diameter'] / 2) < 1e-5 and abs(y) > DATA['derived']['fillet_center'][1]
            is_end = abs(abs(y) - D['overall_length'] / 2) < 1e-5
            if not is_neck and not is_end:
                inside.append(face.index)
    reached, pending = set(), [0]
    while pending:
        index = pending.pop()
        if index not in reached:
            reached.add(index)
            pending.extend(graph[index] - reached)
    assert len(reached) == 3, '金属与塑料之间存在未连接的间隙'
    for label, indices in [('SECTOR_0', side0), ('SECTOR_MINUS_10', side10),
                           ('WCM_OUTER', external), ('INNER_CAVITY', inside),
                           ('MATERIAL_INTERFACES', internal)]:
        assert indices, label + ' 未找到'
        faces = face_array(part, indices)
        part.Set(name=label, faces=faces)
        if label != 'MATERIAL_INTERFACES':
            part.Surface(name=label, side1Faces=faces)
    meridian = [edge.index for edge in part.edges
                if abs(edge.pointOn[0][2]) < 1e-6 and on_wcm(edge.pointOn[0])]
    part.Set(name='WCM_MERIDIAN_0', edges=edge_array(part, meridian))
    axis = part.DatumAxisByPrincipalAxis(principalAxis=YAXIS)
    part.features.changeKey(fromName=axis.name, toName='WCM_REVOLUTION_AXIS_Y')
    sketch_for(model, 'WCM_Meridian', [DATA['wcm_profile']])

    minimum_edge = min(e.getSize(printResults=False) for e in part.edges)
    assert minimum_edge > 1e-4, '存在极短边'
    REPORT['volume_crosschecks'] = []
    for key, actual in [('liner', source_volumes['PA6']), ('boss_positive', source_volumes['Al6061_T6'] / 2)]:
        reference = DATA['checks'][key]['sampled_sector_volume_mm3']
        error = abs(actual - reference) / actual
        REPORT['volume_crosschecks'].append({'region': key, 'abaqus_volume': actual, 'sampled_volume': reference, 'relative_error': error})
        print('VOLUME_CHECK', key, actual, reference, error)
        assert error < 1e-4, '旋转体积与独立截面积分不一致'

    REPORT.update({'status': 'VERIFIED', 'part': part.name,
        'source_volumes_mm3': source_volumes, 'merged_volume_mm3': merged_volume,
        'overlap_volume_loss_mm3': loss, 'relative_volume_loss': relative_loss,
        'cells': len(part.cells), 'connected_cells': len(reached),
        'shared_interface_faces': len(internal), 'sector_angle_range_deg': [min(angles), max(angles)],
        'axial_extent_mm': extent, 'radius_range_mm': [min(radii), max(radii)],
        'minimum_edge_length_mm': minimum_edge,
        'geometry_validity': str(part.geometryValidity),
        'face_sets': {key: len(part.sets[key].faces) for key in ('SECTOR_0', 'SECTOR_MINUS_10', 'WCM_OUTER', 'INNER_CAVITY')},
        'wcm_meridian_edges': len(meridian), 'sections': ['Section_PA6', 'Section_Al6061_T6']})
    assert part.geometryValidity == ON
    cae = os.path.join(RUN, 'zhuning_base_10deg.cae')
    mdb.saveAs(pathName=cae)
    part.writeAcisFile(fileName=os.path.join(RUN, 'zhuning_base_10deg.sat'))
    REPORT['cae_file'] = cae
    write_report()
    with open(os.path.join(ROOT, 'latest_run.json'), 'w', encoding='utf-8') as f:
        json.dump({'run_directory': RUN, 'cae_file': cae}, f, ensure_ascii=False, indent=2)
    print('GEOMETRY_VERIFIED ' + cae)
    print(json.dumps({k: v for k, v in REPORT.items() if k != 'parameters'}, ensure_ascii=True, indent=2))


try:
    main()
except Exception:
    REPORT['status'] = 'FAILED'
    REPORT['error'] = traceback.format_exc()
    write_report()
    print(REPORT['error'])
    raise
