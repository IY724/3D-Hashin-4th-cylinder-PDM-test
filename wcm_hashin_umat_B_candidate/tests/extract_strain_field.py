# -*- coding: utf-8 -*-
"""阶段5-1：从既有fallback ODB只读提取复材part逐积分点应变场。

由 abaqus python 调用：
  abaqus python extract_strain_field.py <odb> <out.json> [frame_time,...]
只读；输出每帧 {time, ips: [[elemLabel, e11,e22,e33,e12,e13,e23], ...]}。
仅导出 TANK-1 实例的 E 场（NLGEOM=NO 下 E 为总机械应变）。
"""
from odbAccess import openOdb
import json
import sys

odb_path = sys.argv[1]
out_path = sys.argv[2]
want_times = [float(x) for x in sys.argv[3].split(',')] if len(sys.argv) > 3 else None

odb = openOdb(odb_path, readOnly=True)
result = {'odb': odb_path, 'frames': []}
try:
    for step in odb.steps.values():
        for frame in step.frames:
            t = frame.frameValue
            if want_times is not None and not any(
                    abs(t-w) < 1e-6 for w in want_times):
                continue
            if 'LE' in frame.fieldOutputs:
                e_field = frame.fieldOutputs['LE']
            elif 'E' in frame.fieldOutputs:
                e_field = frame.fieldOutputs['E']
            else:
                continue
            ips = []
            for value in e_field.values:
                if str(value.instance.name).upper().find('TANK-1') < 0:
                    continue
                try:
                    d = [float(x) for x in value.data]
                except TypeError:
                    d = [float(value.data)]
                ips.append([int(value.elementLabel)] + d[:6])
            result['frames'].append({'time': t, 'ips': ips})
finally:
    odb.close()
with open(out_path, 'w') as handle:
    json.dump(result, handle)
print('frames written:', len(result['frames']),
      'ips per frame:', [len(f['ips']) for f in result['frames']])
