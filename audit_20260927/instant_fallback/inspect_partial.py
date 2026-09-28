from odbAccess import openOdb
from pathlib import Path
path = Path(__file__).resolve().parents[2] / 'test' / 'B_energy_softening' / 'runs' / 'B_20260927_153131_fallback' / 'Rebuild_PDM_B.odb'
odb = openOdb(str(path), readOnly=True)
try:
    for name, step in odb.steps.items():
        print(name, 'timePeriod=', step.timePeriod, 'frames=', len(step.frames))
        print('frame_times=', [f.frameValue for f in step.frames])
finally:
    odb.close()
