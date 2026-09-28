# -*- coding: utf-8 -*-
"""由 abaqus python 调用，只读提取验证单元的积分点数据。

默认保留原行为：每步只取末帧。新增 --all-frames 逐帧输出。
证据保护：输出文件自带来源（odb 路径/大小/修改时间）；目标已存在且
来源不同时拒绝覆盖，防止混装不同运行的证据（规划4.3）。
"""
from odbAccess import openOdb
import json
import os
import sys

args = [a for a in sys.argv[1:] if a not in ('--all-frames', '--overwrite')]
all_frames = '--all-frames' in sys.argv[1:]
overwrite = '--overwrite' in sys.argv[1:]
if len(args) != 2:
    raise SystemExit('用法：extract_odb.py <odb> <output.json> [--all-frames] [--overwrite]')
odb_path, out_path = args[0], args[1]
provenance = {'odb': os.path.abspath(odb_path),
              'odb_bytes': os.path.getsize(odb_path),
              'odb_mtime': os.path.getmtime(odb_path),
              'all_frames': all_frames}
if os.path.exists(out_path) and not overwrite:
    with open(out_path) as _h:
        previous = json.load(_h)
    if not isinstance(previous, dict) or previous.get('provenance') != provenance:
        raise SystemExit('拒绝覆盖来源不同的既有提取结果：' + out_path)

odb = openOdb(odb_path, readOnly=True)
frames = []
for step in odb.steps.values():
    chosen = step.frames if all_frames else (step.frames[-1],)
    for frame in chosen:
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
payload = {'provenance': provenance, 'frames': frames}
with open(out_path, 'w') as handle:
    json.dump(payload, handle)
