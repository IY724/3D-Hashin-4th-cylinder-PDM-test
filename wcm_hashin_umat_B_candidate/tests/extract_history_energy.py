# -*- coding: utf-8 -*-
"""由 abaqus python 调用：提取整模能量历史（ALLSE/ALLIE/ALLPD/ALLAE）。"""
from odbAccess import openOdb
import json
import sys

odb = openOdb(sys.argv[1], readOnly=True)
hist = {}
for step in odb.steps.values():
    for name, region in step.historyRegions.items():
        for key, out in region.historyOutputs.items():
            if key in ('ALLSE', 'ALLIE', 'ALLPD', 'ALLAE'):
                hist.setdefault(key, [float(x[1]) for x in out.data])
odb.close()
with open(sys.argv[2], 'w') as handle:
    json.dump(hist, handle)
print('history keys:', sorted(hist))
