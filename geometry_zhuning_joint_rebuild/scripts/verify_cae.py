# -*- coding: utf-8 -*-
"""保存重开校验：打开 zhuning_joint_rebuilt.cae，逐项核对材料卡与截面绑定。

期望值全部取自 test/A_stiffness_reduction/Job2_PDM_A.inp（与 Job2_PDM_B.inp 一致）。
运行（在 geometry_zhuning_joint_rebuild 目录下）：
    E:\\ABAQUS2025\\Commands\\abaqus.BAT cae noGUI=scripts/verify_cae.py
"""
from abaqus import openMdb
import json
import os
import traceback

CAE_NAME = 'zhuning_joint_rebuilt.cae'
MODEL_NAME = 'Zhuning_Joint_Rebuilt'
PART_NAME = 'Vessel_Base'
EXPECTED = {
    'PA6': {'density': [[1.13e-09]], 'elastic': [[1880.0, 0.4]],
            'plastic': [[47.0, 0.0], [55.0, 0.17]]},
    'al-6061': {'density': [[2.7e-09]], 'elastic': [[70000.0, 0.32]],
                'plastic': [[306.0, 0.0], [340.0, 0.12]]},
    't-700': {'density': [[1.75e-09]],
              'elastic': [[141000.0, 11400.0, 11400.0, 0.28, 0.28, 0.28, 7100.0, 7100.0, 7100.0]],
              'failStrain': [[0.01475, 0.00887, 0.00526, 0.02544, 0.01549]],
              'failStress': [[2080.0, 1250.0, 60.0, 290.0, 110.0, 0.0, 0.0]]},
}


def find_root():
    candidates = []
    if '__file__' in globals():
        candidates.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    candidates.append(os.getcwd())
    for c in candidates:
        if os.path.exists(os.path.join(c, CAE_NAME)):
            return c
    raise RuntimeError('未找到 ' + CAE_NAME)


def region_cell_count(region, part):
    # 重开CAE后 sa.region 是 (集合名, 部件名, ...) 元组，需按集合名回查单元数；
    # 会话内新建时则仍是带 cells 的 Region/Set 对象。
    if hasattr(region, 'cells'):
        return len(region.cells)
    if isinstance(region, (tuple, list)) and region and isinstance(region[0], str):
        set_name = region[0]
        if set_name in part.sets:
            return len(part.sets[set_name].cells)
    return -1


ROOT = find_root()
OUT = os.path.join(ROOT, 'results', 'verify_cae.json')
RESULT = {'status': 'CHECKING', 'cae_file': os.path.join(ROOT, CAE_NAME), 'mismatches': [], 'warnings': []}


def rows(table):
    return [list(r) for r in table]


def check(name, key, actual):
    expected = EXPECTED[name].get(key)
    if expected is None:
        RESULT['warnings'].append('%s.%s 期望值未定义，实际=%s' % (name, key, actual))
        return
    if actual != expected:
        RESULT['mismatches'].append({'material': name, 'item': key, 'expected': expected, 'actual': actual})


def main():
    mdb = openMdb(pathName=RESULT['cae_file'])
    model = mdb.models[MODEL_NAME]
    RESULT['materials_in_model'] = sorted(model.materials.keys())
    for name in EXPECTED:
        mat = model.materials[name]
        check(name, 'density', rows(mat.density.table))
        check(name, 'elastic', rows(mat.elastic.table))
        plastic = getattr(mat, 'plastic', None)
        actual_plastic = rows(plastic.table) if plastic is not None and plastic.table else []
        if 'plastic' in EXPECTED[name]:
            check(name, 'plastic', actual_plastic)
        elif actual_plastic:
            RESULT['mismatches'].append({'material': name, 'item': 'plastic',
                                         'expected': [], 'actual': actual_plastic})
        for key in ('failStrain', 'failStress'):
            if key not in EXPECTED[name]:
                continue
            obj = getattr(mat.elastic, key, None)
            actual = rows(obj.table) if obj is not None and obj.table else []
            check(name, key, actual)

    part = model.parts[PART_NAME]
    RESULT['part'] = PART_NAME
    RESULT['part_cells'] = len(part.cells)
    RESULT['geometry_validity'] = str(part.geometryValidity)
    RESULT['sets'] = sorted(part.sets.keys())
    assignments = []
    for sa in part.sectionAssignments:
        assignments.append({'section': sa.sectionName, 'cells': region_cell_count(sa.region, part)})
    RESULT['section_assignments'] = assignments
    expected_assignments = {'Section_PA6': 1, 'Section_al-6061': 2}
    # 合并体保留源Part映射加脚本显式覆盖（每截面2条），要求每条指派单元数一致即可。
    for a in assignments:
        n = expected_assignments.get(a['section'])
        if n is None or a['cells'] != n:
            RESULT['mismatches'].append({'item': 'section_assignment', 'section': a['section'],
                                         'expected_cells': expected_assignments.get(a['section']),
                                         'actual_cells': a['cells']})
    actual_map = {}
    for a in assignments:
        actual_map[a['section']] = actual_map.get(a['section'], 0) + a['cells']
    RESULT['section_assignment_totals'] = actual_map
    if 'Section_t-700' in actual_map:
        RESULT['mismatches'].append({'item': 'section_assignment', 'section': 'Section_t-700',
                                     'expected_cells': 0, 'actual_cells': actual_map['Section_t-700']})
    RESULT['status'] = 'PASSED' if not RESULT['mismatches'] else 'FAILED'
    RESULT['model_name'] = model.name
    RESULT['assembly_instances'] = sorted(model.rootAssembly.instances.keys())


try:
    main()
except Exception:
    RESULT['status'] = 'ERROR'
    RESULT['error'] = traceback.format_exc()
finally:
    with open(OUT, 'w', encoding='utf-8') as f:
        json.dump(RESULT, f, ensure_ascii=False, indent=2)
    print('VERIFY ' + RESULT['status'] + ' -> ' + OUT)
