# -*- coding: utf-8 -*-
"""Abaqus/CAE 2025 原生建模：10° 基体几何（多 Part 版，不做布尔合并）。

与 scripts/build_abaqus.py 的差异：
- 保留 PA6_Liner 与 Al6061_Bosses 为两个独立 Part（不合并为 Vessel_Base）；
- 材料卡与旧成功案例 G:\\cylinder\\qiping\\workplace\\zhuning\\WCM_Job.inp 一致：
    * PA6     : Elastic 1880/0.4，Plastic 47.@0.0 -> 55.@0.17（两点硬化）
    * Al6061_T6: Elastic 70000/0.32，Plastic 306.@0.0 -> 340.@0.12（两点硬化）
    * T-700   : 原始单层正交各向异性（ENGINEERING CONSTANTS + Fail Stress/Strain），
                非 WCM Bin 等效卡；弹性/失效数值与旧成功案例 T-700 一致，
                密度按工作区记录补 1.75e-9（旧案例 T-700 卡未列密度，不影响静力）；
    密度按本模型单位制 mm-N-MPa-tonne 记 1.13e-9 / 2.7e-9
    （物理值与旧案例 1130 / 2700 kg*mm^-3 相同）；
- 只创建模型与输出目录，不提交分析作业。

运行方式（本脚本位于 scripts\\，数据根目录为其上级）：
    abaqus cae noGUI=scripts/build_multipart_abaqus.py
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

CAE_NAME = 'zhuning_base_10deg_multipart.cae'
MODEL_NAME = 'Zhuning_Base_10deg_Multipart'


def find_root():
    """定位含 profiles.json 的数据根目录（兼容 noGUI 下无 __file__ 的情形）。"""
    candidates = []
    if '__file__' in globals():
        candidates.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    candidates.append(os.getcwd())
    candidates.append(os.path.dirname(os.getcwd()))
    for c in candidates:
        if os.path.exists(os.path.join(c, 'profiles.json')):
            return c
    raise RuntimeError('未找到 profiles.json；请在 geometry_zhuning_10deg 或其 scripts 子目录下运行')


ROOT = find_root()
with open(os.path.join(ROOT, 'profiles.json'), 'r', encoding='utf-8') as f:
    DATA = json.load(f)
P = DATA['parameters']
D = P['paper_dimensions']
RUN = os.path.join(ROOT, 'runs', datetime.datetime.now().strftime('%Y%m%d_%H%M%S') + '_multipart')
os.makedirs(RUN)
os.chdir(RUN)
REPORT = {
    'status': 'BUILDING',
    'run_directory': RUN,
    'model': MODEL_NAME,
    'cae_file': os.path.join(ROOT, CAE_NAME),
    'merged_into_single_part': False,
    'analysis_submitted': False,
    'material_source': 'G:/cylinder/qiping/workplace/zhuning/WCM_Job.inp (旧成功案例)',
    'warnings': [],
}


def write_report():
    with open(os.path.join(RUN, 'multipart_validation.json'), 'w', encoding='utf-8') as f:
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


def main():
    cae_path = os.path.join(ROOT, CAE_NAME)
    if os.path.exists(cae_path):
        raise RuntimeError('目标 CAE 已存在：%s；如需重建请先自行移除该文件。' % cae_path)
    if MODEL_NAME in mdb.models:
        raise RuntimeError('当前会话已存在同名模型；请在新的 CAE 会话运行，避免覆盖已有工作。')
    model = mdb.Model(name=MODEL_NAME)

    # ---- 材料卡：与旧成功案例 WCM_Job.inp 保持一致（含两点硬化） ----
    pa6 = model.Material(name='PA6')
    pa6.Density(table=((1.13e-09,),))          # 与旧案例 1130 (kg/m^3 物理同值)
    pa6.Elastic(table=((1880.0, 0.4),))
    pa6.Plastic(table=((47.0, 0.0), (55.0, 0.17)))
    al = model.Material(name='Al6061_T6')
    al.Density(table=((2.7e-09,),))            # 与旧案例 2700 (kg/m^3 物理同值)
    al.Elastic(table=((70000.0, 0.32),))
    al.Plastic(table=((306.0, 0.0), (340.0, 0.12)))

    # ---- 复合材料原始单层卡：T-700（备用，供后续 WCM/手工铺层引用；无截面指派） ----
    t700 = model.Material(name='T-700')
    t700.Density(table=((1.75e-09,),))         # 工作区记录值；旧案例 T-700 卡未列密度
    t700.Elastic(type=ENGINEERING_CONSTANTS,
                 table=((141000.0, 11400.0, 11400.0, 0.28, 0.28, 0.28,
                         7100.0, 7100.0, 7100.0),))
    # 注意：Abaqus 2025 的 CAE Python API 探测不到 FailStrain/FailStress 方法（2026-09-26 实测，
    # noGUI 下报 AttributeError）；失败卡请在 CAE GUI 手工补录，此处容错避免脚本中断。
    try:
        t700.FailStrain(table=((0.015, 0.009, 0.005, 0.025, 0.015),))
        t700.FailStress(table=((2080.0, 1250.0, 60.0, 290.0, 110.0, 0.0, 110.0),))
    except AttributeError:
        REPORT['warnings'].append(
            'T-700 Fail Stress/Fail Strain 未自动写入（Abaqus 2025 CAE API 不支持），需在 CAE GUI 手工补录')

    model.HomogeneousSolidSection(name='Section_PA6', material='PA6')
    model.HomogeneousSolidSection(name='Section_Al6061_T6', material='Al6061_T6')

    # ---- 两个独立 Part：塑料内衬与金属阀座，不做布尔合并 ----
    liner = make_part(model, 'PA6_Liner', [DATA['liner']], 'PA6', 'MAT_PA6')
    bosses = make_part(model, 'Al6061_Bosses', [DATA['boss_positive'], DATA['boss_negative']],
                       'Al6061_T6', 'MAT_AL6061')
    assert len(liner.cells) == 1, '内衬应为单连通实体'
    assert len(bosses.cells) == 2, '两个铝阀座应各自独立成体'

    liner_volume = volume(liner)
    bosses_volume = volume(bosses)
    assert abs(liner_volume - DATA['checks']['liner']['sampled_sector_volume_mm3']) / liner_volume < 1e-4, \
        '内衬旋转体积与独立截面积分不一致'
    assert abs(bosses_volume / 2.0 - DATA['checks']['boss_positive']['sampled_sector_volume_mm3']) / \
        (bosses_volume / 2.0) < 1e-4, '阀座旋转体积与独立截面积分不一致'

    # ---- 联合几何核对（等价的合并前检查） ----
    coordinates = [v.pointOn[0] for part in (liner, bosses) for v in part.vertices]
    extent_y = (min(q[1] for q in coordinates), max(q[1] for q in coordinates))
    angles = [math.degrees(math.atan2(z, x)) for x, y, z in coordinates]
    radii = [math.hypot(q[0], q[2]) for q in coordinates]
    assert abs(extent_y[1] - extent_y[0] - D['overall_length']) < 1e-6, '联合长度不等于总长 1520'
    assert abs(max(radii) * 2 - D['liner_outer_diameter']) < 1e-6, '最大半径不等于内衬外径'
    assert abs(min(radii) * 2 - D['bore_diameter']) < 1e-6, '最小半径不等于通孔直径'
    assert abs(min(angles) + D['sector_angle_deg']) < 1e-5 and abs(max(angles)) < 1e-5, '扇区角度范围异常'

    # ---- 装配：两个独立实例，不合并 ----
    assembly = model.rootAssembly
    assembly.DatumCsysByDefault(CARTESIAN)
    assembly.Instance(name='PA6_Liner-1', part=liner, dependent=ON)
    assembly.Instance(name='Al6061_Bosses-1', part=bosses, dependent=ON)
    assert len(assembly.instances) == 2, '装配中应保留两个独立实例'

    mdb.saveAs(pathName=cae_path)
    REPORT.update({
        'status': 'VERIFIED',
        'parts': [part.name for part in model.parts.values()],
        'part_cells': {liner.name: len(liner.cells), bosses.name: len(bosses.cells)},
        'materials': {
            'PA6': {'density': 1.13e-09, 'elastic': [1880.0, 0.4],
                    'plastic': [[47.0, 0.0], [55.0, 0.17]]},
            'Al6061_T6': {'density': 2.7e-09, 'elastic': [70000.0, 0.32],
                          'plastic': [[306.0, 0.0], [340.0, 0.12]]},
            'T-700': {'density': 1.75e-09,
                      'elastic_engineering_constants':
                          [141000.0, 11400.0, 11400.0, 0.28, 0.28, 0.28,
                           7100.0, 7100.0, 7100.0],
                      'fail_strain': [0.015, 0.009, 0.005, 0.025, 0.015],
                      'fail_stress': [2080.0, 1250.0, 60.0, 290.0, 110.0, 0.0, 110.0]},
        },
        'sections': ['Section_PA6 -> MAT_PA6', 'Section_Al6061_T6 -> MAT_AL6061'],
        'instances': ['PA6_Liner-1', 'Al6061_Bosses-1'],
        'volumes_mm3': {'PA6_Liner': liner_volume, 'Al6061_Bosses_total': bosses_volume,
                        'Al6061_Bosses_per_boss': bosses_volume / 2.0},
        'axial_extent_mm': extent_y,
        'radius_range_mm': [min(radii), max(radii)],
        'geometry_validity': {'PA6_Liner': str(liner.geometryValidity),
                              'Al6061_Bosses': str(bosses.geometryValidity)},
    })
    write_report()
    print('MULTIPART_CAE_VERIFIED ' + cae_path)
    print(json.dumps({k: v for k, v in REPORT.items() if k != 'parameters'}, ensure_ascii=True, indent=2))


try:
    main()
except Exception:
    REPORT['status'] = 'FAILED'
    REPORT['error'] = traceback.format_exc()
    write_report()
    print(REPORT['error'])
    raise
