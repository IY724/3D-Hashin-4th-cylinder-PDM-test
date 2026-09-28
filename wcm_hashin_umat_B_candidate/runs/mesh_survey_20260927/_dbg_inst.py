from odbAccess import openOdb
odb = openOdb(r'..\test\B_energy_softening\runs\B_20260927_153131_fallback\Rebuild_PDM_B.odb', True)
step = odb.steps.values()[0]
fr = step.frames[8]
print('t=', fr.frameValue)
print('instances:', [inst.name for inst in odb.rootAssembly.instances.values()])
print('fields:', [k for k in fr.fieldOutputs.keys()])
e = fr.fieldOutputs['E']
vs = e.values
print('E values:', len(vs), 'first instance:', vs[0].instance.instanceName)
locs = set(str(v.instance.instanceName) for v in vs[:5000])
print('E instances sample:', locs)
print('E positions:', set(str(v.position) for v in vs[:2000]))
odb.close()
