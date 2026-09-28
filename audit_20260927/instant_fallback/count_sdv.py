from odbAccess import openOdb
from pathlib import Path
import json

root = Path(__file__).resolve().parents[2]
path = root / 'test' / 'B_energy_softening' / 'runs' / 'B_20260927_153131_fallback' / 'Rebuild_PDM_B.odb'
odb = openOdb(str(path), readOnly=True)
try:
    step = odb.steps['Step-1']
    results = []
    for frame in step.frames:
        if not all('SDV%d' % i in frame.fieldOutputs for i in (5, 6, 7, 8)):
            continue
        row = {'step_time': frame.frameValue}
        for i in (5, 6, 7, 8):
            values = frame.fieldOutputs['SDV%d' % i].values
            scalars = [float(v.data) for v in values]
            row['sdv%d_at_one' % i] = sum(v >= 1.0-1e-6 for v in scalars)
            row['sdv%d_max' % i] = max(scalars) if scalars else None
            row['points'] = len(scalars)
        results.append(row)
finally:
    odb.close()
target = Path(__file__).with_name('sdv_counts.json')
target.write_text(json.dumps(results, indent=2))
print(str(target))
