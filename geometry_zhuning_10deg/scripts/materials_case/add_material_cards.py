# -*- coding: utf-8 -*-
"""按用户要求原位添加旧案例材料卡；不指派截面、不创建铺层、不另存CAE。"""
from abaqus import mdb, openMdb
from abaqusConstants import *
import caeModules
from builtins import all, min, max
import hashlib
import json
import os
import traceback

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) if '__file__' in globals() else os.getcwd()  # 数据根目录为 scripts\materials_case 的上两级
with open(os.path.join(ROOT, 'case_materials_and_layup.json'), 'r', encoding='utf-8') as f:
    SOURCE = json.load(f)
with open(os.path.join(ROOT, 'latest_run.json'), 'r', encoding='utf-8') as f:
    CAE = json.load(f)['cae_file']
MODEL = 'Zhuning_Base_10deg'
STATUS = {'status': 'RUNNING', 'cae_file': CAE, 'save_in_place': True,
          'new_section_assignments': False, 'new_layups': False, 'analysis_submitted': False}


def file_hash(path):
    with open(path, 'rb') as f:
        return hashlib.sha256(f.read()).hexdigest()


def write_status():
    with open(os.path.join(ROOT, 'material_update_validation.json'), 'w', encoding='utf-8') as f:
        json.dump(STATUS, f, ensure_ascii=False, indent=2)


def material_snapshot(mat):
    result = {}
    for attribute in ('density', 'elastic', 'plastic'):
        if hasattr(mat, attribute):
            value = getattr(mat, attribute)
            result[attribute] = [list(row) for row in value.table]
    if hasattr(mat, 'elastic'):
        result['elastic_type'] = str(mat.elastic.type)
        for attribute in ('failStress', 'failStrain'):
            if hasattr(mat.elastic, attribute):
                result[attribute] = [list(row) for row in getattr(mat.elastic, attribute).table]
    return result


def configuration_snapshot(model):
    parts = {}
    for name, part in model.parts.items():
        parts[name] = {'counts': [len(part.cells), len(part.faces), len(part.edges), len(part.vertices),
                                 len(part.nodes), len(part.elements)],
                       'vertices': [list(v.pointOn[0]) for v in part.vertices],
                       'sets': sorted(part.sets.keys()), 'surfaces': sorted(part.surfaces.keys()),
                       'features': sorted(part.features.keys()),
                       'layups': sorted(part.compositeLayups.keys()),
                       'assignments': [(sa.region[0], sa.sectionName) for sa in part.sectionAssignments]}
    return {'parts': parts, 'sections': {n: s.material for n, s in model.sections.items()},
            'steps': sorted(model.steps.keys()), 'loads': sorted(model.loads.keys()),
            'boundary_conditions': sorted(model.boundaryConditions.keys()),
            'constraints': sorted(model.constraints.keys()),
            'assembly_features': sorted(model.rootAssembly.features.keys()),
            'instances': sorted(model.rootAssembly.instances.keys()),
            'jobs': sorted(mdb.jobs.keys()), 'models': sorted(mdb.models.keys())}


def desired_properties(card):
    result = {}
    for prop in card['properties']:
        key, values = prop['keyword'], prop['values']
        if key == 'elastic':
            engineering = prop['options'].get('type') == 'ENGINEERING CONSTANTS'
            result['elastic_type'] = 'ENGINEERING_CONSTANTS' if engineering else 'ISOTROPIC'
            result['elastic'] = [[v for row in values for v in row]] if engineering else values
        else:
            attribute = {'density': 'density', 'plastic': 'plastic',
                         'fail strain': 'failStrain', 'fail stress': 'failStress'}[key]
            result[attribute] = values
    return result


try:
    if os.path.exists(os.path.splitext(CAE)[0] + '.lck'):
        raise RuntimeError('CAE正在被其他会话使用，未修改文件。请先保存并关闭该数据库。')
    STATUS['cae_sha256_before'] = file_hash(CAE)
    for source in SOURCE['sources']:
        assert file_hash(source['path']) == source['sha256'], '来源INP已变化，请先重新提取'
    openMdb(pathName=CAE)
    model = mdb.models[MODEL]
    before = configuration_snapshot(model)
    existing = {name: material_snapshot(mat) for name, mat in model.materials.items()}
    targets = SOURCE['materials']
    expected = {name: desired_properties(card) for name, card in targets.items()}
    created = []
    for name, values in expected.items():
        if name in model.materials:
            assert material_snapshot(model.materials[name]) == values, '已有同名材料参数不同，停止以避免覆盖'
            continue
        mat = model.Material(name=name)
        if 'density' in values:
            mat.Density(table=tuple(tuple(row) for row in values['density']))
        mat.Elastic(type=ENGINEERING_CONSTANTS if values['elastic_type'] == 'ENGINEERING_CONSTANTS' else ISOTROPIC,
                    table=tuple(tuple(row) for row in values['elastic']))
        if 'plastic' in values:
            mat.Plastic(table=tuple(tuple(row) for row in values['plastic']))
        if 'failStrain' in values:
            mat.elastic.FailStrain(table=tuple(tuple(row) for row in values['failStrain']))
        if 'failStress' in values:
            mat.elastic.FailStress(table=tuple(tuple(row) for row in values['failStress']))
        assert material_snapshot(mat) == values, name + '参数不一致'
        created.append(name)
    assert configuration_snapshot(model) == before, '材料卡之外的配置发生变化'
    assert all(material_snapshot(model.materials[n]) == v for n, v in existing.items()), '原有材料被改变'
    mdb.save()
    mdb.close()
    openMdb(pathName=CAE)
    model = mdb.models[MODEL]
    assert configuration_snapshot(model) == before, '保存重开后配置不一致'
    for name, values in expected.items():
        assert material_snapshot(model.materials[name]) == values, '材料卡未正确持久化'
    assert all(material_snapshot(model.materials[n]) == v for n, v in existing.items())
    STATUS.update({'status': 'VERIFIED', 'created_materials': created,
                   'materials': {n: material_snapshot(model.materials[n]) for n in expected},
                   'preserved_existing_materials': sorted(existing),
                   'geometry_and_configuration_unchanged': True,
                   'round_trip_reopen_verified': True,
                   'source_inputs_unchanged': all(file_hash(q['path']) == q['sha256'] for q in SOURCE['sources'])})
    mdb.close()
    STATUS['cae_sha256_after'] = file_hash(CAE)
    write_status()
    print('MATERIAL_CARDS_VERIFIED_IN_PLACE: ' + CAE)
except Exception:
    STATUS['status'] = 'FAILED'
    STATUS['error'] = traceback.format_exc()
    write_status()
    raise
