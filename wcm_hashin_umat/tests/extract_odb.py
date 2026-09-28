# -*- coding: utf-8 -*-
"""由 abaqus python 调用，只读提取验证单元的积分点数据。"""
from odbAccess import openOdb
import json
import sys

odb = openOdb(sys.argv[1], readOnly=True)
frames = []
for step in odb.steps.values():
    for frame in (step.frames[-1],):
        row = {'time': frame.frameValue, 'fields': {}}
        for name in frame.fieldOutputs.keys():
            if name not in ('S', 'E', 'LE') and not name.startswith('SDV'):
                continue
            entries = []
            for value in frame.fieldOutputs[name].values:
                try:
                    data = [float(x) for x in value.data]
                except TypeError:
                    data = [float(value.data)]
                entries.append({'element': value.elementLabel,
                                'point': getattr(value, 'integrationPoint', 0),
                                'data': list(data)})
            row['fields'][name] = entries
        frames.append(row)
odb.close()
with open(sys.argv[2], 'w') as handle:
    json.dump(frames, handle)
