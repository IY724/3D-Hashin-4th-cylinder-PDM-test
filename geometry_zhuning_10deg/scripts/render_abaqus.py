# -*- coding: utf-8 -*-
"""直接渲染已验证的 CAE 实体，不调用任何 AI 图像生成服务。

推荐命令：abaqus cae noGUI=render_abaqus.py
保存截面、端部放大及10°实体透视图。只读取已生成 CAE。
"""
from abaqus import mdb, session, openMdb
from abaqusConstants import *
import caeModules
import os
import json
import traceback

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) if '__file__' in globals() else os.getcwd()  # 数据根目录为 scripts\ 的上级
with open(os.path.join(ROOT, 'latest_run.json'), 'r', encoding='utf-8') as f:
    latest = json.load(f)
OUT = latest['run_directory']
STATUS = {'source_cae': latest['cae_file'], 'renderer': 'Abaqus/CAE native viewport',
          'ai_generated': False, 'images': [], 'status': 'RUNNING'}


def save_status():
    with open(os.path.join(OUT, 'render_status.json'), 'w', encoding='utf-8') as f:
        json.dump(STATUS, f, ensure_ascii=False, indent=2)


try:
    openMdb(pathName=latest['cae_file'])
    part = mdb.models['Zhuning_Base_10deg'].parts['Vessel_Base']
    part.setValues(geometryRefinement=EXTRA_FINE)
    vp = session.viewports[session.currentViewportName]
    vp.setValues(displayedObject=part)
    vp.partDisplay.setValues(renderStyle=SHADED)
    vp.partDisplay.geometryOptions.setValues(datumAxes=OFF, datumPlanes=OFF,
                                             datumPoints=OFF, datumCoordSystems=OFF)
    session.graphicsOptions.setValues(backgroundStyle=SOLID, backgroundColor='#FFFFFF')
    vp.viewportAnnotationOptions.setValues(triad=ON, legend=OFF, title=OFF, state=OFF)
    vp.enableMultipleColors()
    vp.setColor(initialColor='#A8A8A8')
    cmap = vp.colorMappings['Material']
    cmap.updateOverrides(overrides={
        'PA6': (True, '#398FBD', 'Default', '#398FBD'),
        'Al6061_T6': (True, '#D9A14C', 'Default', '#D9A14C')})
    vp.setColor(colorMapping=cmap)
    vp.disableMultipleColors()
    session.printOptions.setValues(vpDecorations=OFF, vpBackground=ON)
    session.pngOptions.setValues(imageSize=(2400, 1400))

    def export(name):
        vp.forceRefresh()
        path = os.path.join(OUT, name)
        session.printToFile(fileName=path, format=PNG, canvasObjects=(vp,))
        STATUS['images'].append(path + '.png')
        save_status()

    # 从正Z方向观察实体的Z=0子午侧面，即真实几何截面。
    vp.view.setValues(projection=PARALLEL, cameraPosition=(70.0, 0.0, 2500.0),
                      cameraTarget=(70.0, 0.0, 0.0), cameraUpVector=(1.0, 0.0, 0.0))
    vp.view.fitView()
    export('abaqus_section_full')
    vp.view.setValues(cameraPosition=(75.0, 680.0, 2500.0),
                      cameraTarget=(75.0, 680.0, 0.0), cameraUpVector=(1.0, 0.0, 0.0),
                      width=260.0, height=170.0)
    export('abaqus_section_joint')
    vp.view.setValues(cameraPosition=(850.0, 850.0, 2500.0),
                      cameraTarget=(70.0, 0.0, -8.0), cameraUpVector=(1.0, 0.0, 0.0))
    vp.view.fitView()
    export('abaqus_sector_3d')
    STATUS['status'] = 'COMPLETE'
    save_status()
except Exception:
    STATUS['status'] = 'FAILED'
    STATUS['error'] = traceback.format_exc()
    save_status()
    raise
