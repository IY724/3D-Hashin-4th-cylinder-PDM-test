# -*- coding: utf-8 -*-
"""Update zhuning_base_10deg_multipart.cae in place: add original T-700 ply card.

- Writes the ORIGINAL single-ply t-700 card (ENGINEERING CONSTANTS + Fail
  Stress/Strain), NOT the WCM Bin smeared cards;
- Values match the legacy successful case WCM_Job.inp T-700; density 1.75e-9
  taken from repository records (legacy card had no density; no static effect);
- Idempotent: skips creation if material already present;
- All log lines are ASCII to avoid GBK console encoding failures in noGUI.

Run (from geometry_zhuning_10deg):
    abaqus cae noGUI=scripts/add_t700_to_multipart.py
"""
import io
import os
import traceback

from abaqus import openMdb
from abaqusConstants import ENGINEERING_CONSTANTS

ROOT = r'f:\abaqus_temp\3D-Hashin-4th-cylinder-PDM\geometry_zhuning_10deg'
CAE = os.path.join(ROOT, 'zhuning_base_10deg_multipart.cae')
OUT = os.path.join(ROOT, 'scripts', '_add_t700_to_multipart.out')

lines = []


def emit(text):
    lines.append(str(text))
    print(text)


def dump_materials(tag):
    m = mdb.models['Zhuning_Base_10deg_Multipart']
    for mn in list(m.materials.keys()):
        mat = m.materials[mn]
        plastic = mat.plastic.table if hasattr(mat, 'plastic') else None
        fail_stress = mat.failStress.table if hasattr(mat, 'failStress') else None
        fail_strain = mat.failStrain.table if hasattr(mat, 'failStrain') else None
        emit('[%s] MATERIAL %s: density=%s elasticType=%s elastic=%s plastic=%s failStress=%s failStrain=%s' %
             (tag, mn, mat.density.table, getattr(mat, 'elasticType', None),
              mat.elastic.table, plastic, fail_stress, fail_strain))


status = 'UNKNOWN'
try:
    openMdb(pathName=CAE)
    emit('STEP openMdb OK')
    m = mdb.models['Zhuning_Base_10deg_Multipart']
    emit('STEP model OK: %s' % m.name)

    if 'T-700' in list(m.materials.keys()):
        emit('STEP T-700 already present; skip creation')
        status = 'ALREADY_PRESENT'
    else:
        t700 = m.Material(name='T-700')
        emit('STEP material created')
        t700.Density(table=((1.75e-09,),))
        t700.Elastic(type=ENGINEERING_CONSTANTS,
                     table=((141000.0, 11400.0, 11400.0, 0.28, 0.28, 0.28,
                             7100.0, 7100.0, 7100.0),))
        t700.FailStrain(table=((0.015, 0.009, 0.005, 0.025, 0.015),))
        t700.FailStress(table=((2080.0, 1250.0, 60.0, 290.0, 110.0, 0.0, 110.0),))
        emit('STEP cards set')
        try:
            mdb.save()
            emit('STEP mdb.save OK')
        except Exception:
            emit('STEP mdb.save failed:\n' + traceback.format_exc())
            mdb.saveAs(pathName=CAE)
            emit('STEP mdb.saveAs fallback OK')
        status = 'WRITTEN'

    dump_materials('memory')
except Exception:
    status = 'FAILED'
    emit('EXCEPTION:\n' + traceback.format_exc())

emit('STATUS=%s' % status)
with io.open(OUT, 'w', encoding='utf-8') as f:
    f.write(u'\n'.join(lines) + u'\n')
