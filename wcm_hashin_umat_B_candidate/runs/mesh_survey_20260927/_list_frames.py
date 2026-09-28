from odbAccess import openOdb
odb = openOdb(r'..\test\B_energy_softening\runs\B_20260927_153131_fallback\Rebuild_PDM_B.odb', True)
for step in odb.steps.values():
    ts = [round(f.frameValue, 6) for f in step.frames]
    print('step', step.name, 'nframes', len(ts))
    print(ts)
odb.close()
