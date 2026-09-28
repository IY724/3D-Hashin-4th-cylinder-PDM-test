# -*- coding: utf-8 -*-
"""Abaqus/CAE 2025 原生建模：朱宁宽肩接头 paper-fit v3 新尺寸，单 Part 交付。

- 截面数据来自 archive/paper_fit_preview/preview_profiles.json（prepare_preview.py 产物），
  尺寸约束见 parameters.json（zhuning-wide-shoulder-paper-fit-v3）；
- 材料卡取自 test/A_stiffness_reduction/Job2_PDM_A.inp 与 test/B_energy_softening/Job2_PDM_B.inp
  （两套 INP 逐项一致）：PA6 / al-6061 / t-700，含 t-700 的 Fail Strain / Fail Stress；
- 两个源 Part（PA6_Liner、Al6061_Bosses）布尔合并为单 Part Vessel_Base（3 cell）：
  PA6 内衬 1 cell、Al6061 阀座 2 cell，合并后显式保留两种材料分区；
- 只建几何与材料分区，不建铺层、网格、分析步、载荷与求解。

运行（在 geometry_zhuning_joint_rebuild 目录下）：
    E:\\ABAQUS2025\\Commands\\abaqus.BAT cae noGUI=scripts/build_single_part_abaqus.py
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

import numpy as np

CAE_NAME = 'zhuning_joint_rebuilt.cae'
MODEL_NAME = 'Zhuning_Joint_Rebuilt'
MATERIAL_SOURCE = ('test/A_stiffness_reduction/Job2_PDM_A.inp 与 test/B_energy_softening/Job2_PDM_B.inp'
                   '（两套一致：PA6 @275987、al-6061 @288440/288571、t-700 @288448/288579）')


def find_root():
    """定位含 archive/paper_fit_preview/preview_profiles.json 的数据根目录。"""
    candidates = []
    if '__file__' in globals():
        candidates.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    candidates.append(os.getcwd())
    candidates.append(os.path.dirname(os.getcwd()))
    for c in candidates:
        if os.path.exists(os.path.join(c, 'archive', 'paper_fit_preview', 'preview_profiles.json')):
            return c
    raise RuntimeError('未找到 preview_profiles.json；请在 geometry_zhuning_joint_rebuild 下运行')


ROOT = find_root()
with open(os.path.join(ROOT, 'archive', 'paper_fit_preview', 'preview_profiles.json'), 'r', encoding='utf-8') as f:
    DATA = json.load(f)
with open(os.path.join(ROOT, 'parameters.json'), 'r', encoding='utf-8') as f:
    P = json.load(f)
D = P['paper']
SECTOR = D['sector_angle']
OVERALL = D['overall_length']
OUTER_D = D['outer_diameter']
BORE_D = D['bore_diameter']
RUN = os.path.join(ROOT, 'archive', 'build_runs',
                   datetime.datetime.now().strftime('%Y%m%d_%H%M%S') + '_singlepart')
os.makedirs(RUN)
os.chdir(RUN)
REPORT = {
    'status': 'BUILDING',
    'run_directory': RUN,
    'model': MODEL_NAME,
    'geometry_source': 'archive/paper_fit_preview/preview_profiles.json (paper-fit v3)',
    'material_source': MATERIAL_SOURCE,
    'merged_into_single_part': True,
    'analysis_submitted': False,
    'layup_created': False,
    'warnings': [],
}


def write_report():
    with open(os.path.join(RUN, 'build_validation.json'), 'w', encoding='utf-8') as f:
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


def make_part(model, name, loops, section, set_name):
    sk = sketch_for(model, 'Sketch_' + name, loops)
    part = model.Part(name=name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
    part.BaseSolidRevolve(sketch=sk, angle=SECTOR, flipRevolveDirection=ON)
    region = part.Set(name=set_name, cells=part.cells)
    part.SectionAssignment(region=region, sectionName=section)
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


# 轮廓段按 prepare_preview.py 的构造顺序标注：OUT=外表面 IN=内腔 END=端面 INT=材料界面。
LINER_LABELS = ['OUT', 'OUT', 'INT', 'INT', 'INT', 'INT',
                'IN', 'IN', 'IN', 'IN', 'IN', 'IN',
                'IN', 'IN', 'IN', 'IN', 'IN',
                'INT', 'INT', 'INT', 'INT', 'OUT']
BOSS_LABELS = ['OUT', 'OUT', 'OUT', 'END',
               'IN', 'IN', 'IN', 'IN', 'IN',
               'INT', 'INT', 'INT', 'INT']


def label_points(entries):
    pts = {'IN': [], 'OUT': [], 'END': []}
    for loop, label in zip([e[0] for e in entries], [e[1] for e in entries]):
        if label in pts:
            # 轮廓点坐标即 (半径, 轴向)，与回转后面点的 (hypot(x,z), y) 同一量纲。
            # line 段只有两个端点，须加密采样至间距不大于 0.2 mm，保证面心点能命中自身标签。
            coords = loop['points']
            for (x0, y0), (x1, y1) in zip(coords[:-1], coords[1:]):
                n = max(1, int(math.ceil(math.hypot(x1 - x0, y1 - y0) / 0.2)))
                for k in range(n):
                    t = k / float(n)
                    pts[label].append((x0 + (x1 - x0) * t, y0 + (y1 - y0) * t))
            pts[label].append(tuple(coords[-1]))
    return {k: np.array(v) for k, v in pts.items()}


assert len(DATA['liner']) == len(LINER_LABELS), '内衬轮廓段数与标签数不一致'
assert len(DATA['boss_positive']) == len(BOSS_LABELS), '阀座轮廓段数与标签数不一致'
assert len(DATA['boss_negative']) == len(BOSS_LABELS), '镜像阀座轮廓段数与标签数不一致'
LABEL_PTS = label_points(
    [(loop, lab) for loop, lab in zip(DATA['liner'], LINER_LABELS)] +
    [(loop, lab) for loop, lab in zip(DATA['boss_positive'], BOSS_LABELS)] +
    [(loop, lab) for loop, lab in zip(DATA['boss_negative'], BOSS_LABELS)])


def classify(q):
    scores = {}
    for label, pts in LABEL_PTS.items():
        scores[label] = float(np.min(np.hypot(pts[:, 0] - q[0], pts[:, 1] - q[1])))
    best = min(scores, key=scores.get)
    # 0.5 mm 容差：远大于采样间距，又小于空腔中心到轮廓的距离（约 4 mm）。
    return best if scores[best] < 0.5 else 'UNCLASSIFIED'


def add_materials(model):
    # 三张卡数值逐项取自 test A/B INP（单位 mm-N-MPa-tonne）。
    pa6 = model.Material(name='PA6')
    pa6.Density(table=((1.13e-09,),))
    pa6.Elastic(table=((1880.0, 0.4),))
    pa6.Plastic(table=((47.0, 0.0), (55.0, 0.17)))
    al = model.Material(name='al-6061')
    al.Density(table=((2.7e-09,),))
    al.Elastic(table=((70000.0, 0.32),))
    al.Plastic(table=((306.0, 0.0), (340.0, 0.12)))
    t7 = model.Material(name='t-700')
    t7.Density(table=((1.75e-09,),))
    t7.Elastic(type=ENGINEERING_CONSTANTS,
               table=((141000.0, 11400.0, 11400.0, 0.28, 0.28, 0.28,
                       7100.0, 7100.0, 7100.0),))
    # FailStrain/FailStress 挂在 elastic 子对象上（2026-09-26 实测），根对象调用会 AttributeError。
    fail_strain = ((0.01475, 0.00887, 0.00526, 0.02544, 0.01549),)
    fail_stress = ((2080.0, 1250.0, 60.0, 290.0, 110.0, 0.0, 0.0),)
    try:
        t7.elastic.FailStrain(table=fail_strain)
        t7.elastic.FailStress(table=fail_stress)
    except Exception as exc:
        REPORT['warnings'].append(
            't-700 Fail Strain/Stress 未自动写入（Abaqus 2025 CAE API 限制：%s）；'
            '数值已记录于 materials 明细，可在 CAE GUI 手工补录' % exc)
    model.HomogeneousSolidSection(name='Section_PA6', material='PA6')
    model.HomogeneousSolidSection(name='Section_al-6061', material='al-6061')
    model.HomogeneousSolidSection(name='Section_t-700', material='t-700')


def main():
    cae_path = os.path.join(ROOT, CAE_NAME)
    if os.path.exists(cae_path):
        raise RuntimeError('目标 CAE 已存在：%s；如需重建请先自行移除该文件。' % cae_path)
    if MODEL_NAME in mdb.models:
        raise RuntimeError('当前会话已存在同名模型；请在新的 CAE 会话运行，避免覆盖已有工作。')
    model = mdb.Model(name=MODEL_NAME)
    add_materials(model)

    liner = make_part(model, 'PA6_Liner', [DATA['liner']], 'Section_PA6', 'MAT_PA6')
    bosses = make_part(model, 'Al6061_Bosses', [DATA['boss_positive'], DATA['boss_negative']],
                       'Section_al-6061', 'MAT_AL6061')
    assert len(liner.cells) == 1, '内衬应为单连通实体'
    assert len(bosses.cells) == 2, '两个铝阀座应各自独立成体'
    liner_volume = volume(liner)
    bosses_volume = volume(bosses)

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
    part.SectionAssignment(region=part.sets['MAT_AL6061'], sectionName='Section_al-6061')

    # ---- 体积与外形核对 ----
    merged_volume = volume(part)
    loss = liner_volume + bosses_volume - merged_volume
    relative_loss = abs(loss) / (liner_volume + bosses_volume)
    assert relative_loss < 1e-8, '合并发生体积损失，可能存在穿模'
    assert abs(liner_volume - DATA['checks']['liner']['volume_mm3']) / liner_volume < 1e-4, \
        '内衬旋转体积与独立截面积分不一致'
    assert abs(bosses_volume / 2.0 - DATA['checks']['boss']['volume_mm3']) / (bosses_volume / 2.0) < 1e-4, \
        '阀座旋转体积与独立截面积分不一致'

    coordinates = [v.pointOn[0] for v in part.vertices]
    extent = [min(q[1] for q in coordinates), max(q[1] for q in coordinates)]
    angles = [math.degrees(math.atan2(q[2], q[0])) for q in coordinates]
    radii = [math.hypot(q[0], q[2]) for q in coordinates]
    assert abs(extent[1] - extent[0] - OVERALL) < 1e-6, '轴向总长不等于 1520'
    assert abs(max(radii) * 2 - OUTER_D) < 1e-6, '最大半径不等于外径 280'
    assert abs(min(radii) * 2 - BORE_D) < 1e-6, '最小半径不等于通孔直径 32'
    assert abs(min(angles) + SECTOR) < 1e-5 and abs(max(angles)) < 1e-5, '扇区角度范围异常'

    # ---- 面分类：材料界面 / 扇区面 / 外表面 / 内腔 / 端面 ----
    graph = {c.index: set() for c in part.cells}
    internal, side0, side10 = [], [], []
    labelled = {'OUT': [], 'IN': [], 'END': [], 'UNCLASSIFIED': []}
    for face in part.faces:
        adjacent = tuple(face.getCells())
        if len(adjacent) == 2:
            internal.append(face.index)
            graph[adjacent[0]].add(adjacent[1])
            graph[adjacent[1]].add(adjacent[0])
            continue
        xyz = face.pointOn[0]
        normal = face.getNormal(point=xyz)
        theta = math.atan2(xyz[2], xyz[0])
        circumferential = (-math.sin(theta), 0.0, math.cos(theta))
        if abs(sum(a * b for a, b in zip(normal, circumferential))) > .99999:
            (side0 if abs(xyz[2]) < 1e-6 else side10).append(face.index)
            continue
        label = classify((math.hypot(xyz[0], xyz[2]), xyz[1]))
        labelled[label].append(face.index)
    reached, pending = set(), [0]
    while pending:
        index = pending.pop()
        if index not in reached:
            reached.add(index)
            pending.extend(graph[index] - reached)
    assert len(reached) == 3, '金属与塑料之间存在未连接的间隙'
    for label, indices in [('SECTOR_0', side0), ('SECTOR_MINUS_10', side10),
                           ('WCM_OUTER', labelled['OUT']), ('INNER_CAVITY', labelled['IN']),
                           ('END_FACES', labelled['END']), ('MATERIAL_INTERFACES', internal)]:
        assert indices, label + ' 未找到'
        faces = face_array(part, indices)
        part.Set(name=label, faces=faces)
        if label != 'MATERIAL_INTERFACES':
            part.Surface(name=label, side1Faces=faces)

    meridian = [edge.index for edge in part.edges
                if abs(edge.pointOn[0][2]) < 1e-6
                and classify((math.hypot(edge.pointOn[0][0], edge.pointOn[0][2]), edge.pointOn[0][1])) == 'OUT']
    assert meridian, '未找到外表面子午线'
    part.Set(name='WCM_MERIDIAN_0', edges=edge_array(part, meridian))
    axis = part.DatumAxisByPrincipalAxis(principalAxis=YAXIS)
    part.features.changeKey(fromName=axis.name, toName='WCM_REVOLUTION_AXIS_Y')

    minimum_edge = min(e.getSize(printResults=False) for e in part.edges)
    assert minimum_edge > 1e-4, '存在极短边'
    assert part.geometryValidity == ON

    REPORT.update({
        'status': 'VERIFIED',
        'part': part.name,
        'part_cells': len(part.cells),
        'sections': ['Section_PA6 -> MAT_PA6 (1 cell)', 'Section_al-6061 -> MAT_AL6061 (2 cells)',
                     'Section_t-700（未指派，备用单层卡）'],
        'materials': {
            'PA6': {'density': 1.13e-09, 'elastic': [1880.0, 0.4],
                    'plastic': [[47.0, 0.0], [55.0, 0.17]]},
            'al-6061': {'density': 2.7e-09, 'elastic': [70000.0, 0.32],
                        'plastic': [[306.0, 0.0], [340.0, 0.12]]},
            't-700': {'density': 1.75e-09,
                      'elastic_engineering_constants': [141000.0, 11400.0, 11400.0, 0.28, 0.28, 0.28,
                                                        7100.0, 7100.0, 7100.0],
                      'fail_strain': [0.01475, 0.00887, 0.00526, 0.02544, 0.01549],
                      'fail_stress': [2080.0, 1250.0, 60.0, 290.0, 110.0, 0.0, 0.0]},
        },
        'volumes_mm3': {'PA6_Liner': liner_volume, 'Al6061_Bosses_total': bosses_volume,
                        'merged_single_part': merged_volume},
        'reference_volumes_mm3': DATA['checks'],
        'overlap_volume_loss_mm3': loss, 'relative_volume_loss': relative_loss,
        'connected_cells': len(reached),
        'shared_interface_faces': len(internal),
        'face_sets': {key: len(part.sets[key].faces) for key in
                      ('SECTOR_0', 'SECTOR_MINUS_10', 'WCM_OUTER', 'INNER_CAVITY', 'END_FACES')},
        'unclassified_faces': labelled['UNCLASSIFIED'],
        'axial_extent_mm': extent,
        'radius_range_mm': [min(radii), max(radii)],
        'sector_angle_range_deg': [min(angles), max(angles)],
        'minimum_edge_length_mm': minimum_edge,
        'geometry_validity': str(part.geometryValidity),
        'wcm_meridian_edges': len(meridian),
    })
    cae = os.path.join(ROOT, CAE_NAME)
    mdb.saveAs(pathName=cae)
    part.writeAcisFile(fileName=os.path.join(RUN, 'zhuning_joint_rebuilt.sat'))
    REPORT['cae_file'] = cae
    write_report()
    status_path = os.path.join(ROOT, 'results', 'paper_fit_preview_status.json')
    if os.path.exists(status_path):
        with open(status_path, 'r', encoding='utf-8') as f:
            status = json.load(f)
        status.update({'cae_saved': True, 'cae_file': cae,
                       'build_report': os.path.join(RUN, 'build_validation.json')})
        with open(status_path, 'w', encoding='utf-8') as f:
            json.dump(status, f, ensure_ascii=False, indent=2)
    print('SINGLEPART_CAE_VERIFIED ' + cae)
    print(json.dumps({k: v for k, v in REPORT.items() if k != 'materials'}, ensure_ascii=True, indent=2))


try:
    main()
except Exception:
    REPORT['status'] = 'FAILED'
    REPORT['error'] = traceback.format_exc()
    write_report()
    print(REPORT['error'])
    raise
