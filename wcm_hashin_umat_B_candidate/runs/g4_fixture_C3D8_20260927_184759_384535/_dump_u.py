from odbAccess import openOdb
odb = openOdb('fixture.odb', True)
step = odb.steps.values()[0]
fr = step.frames[-3]
print('frame time', fr.frameValue)
fu = fr.fieldOutputs['U']
for v in fu.values:
    print('node', v.nodeLabel, list(v.data))
odb.close()
